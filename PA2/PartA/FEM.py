import sys
import os
import numpy as np
import meshio

# Sim Constants
c0 = 299792458.0
mu0 = 4.0e-7 * np.pi
eps0 = 1.0 / (mu0 * c0 * c0)
pecTag = 1000000


def loadMesh(meshFile):
    mesh = meshio.read(meshFile)
    nodeXY = np.asarray(mesh.points)[:, :2].astype(float)
    tri = np.asarray(mesh.cells_dict["triangle"], dtype=np.int64)
    tag = pecTag
    for name, val in (getattr(mesh, "field_data", {}) or {}).items():
        if str(name).lower() == "pec":
            tag = int(np.asarray(val).ravel()[0])
    lineElems = np.concatenate(
        [np.asarray(cb.data, dtype=np.int64) for cb in mesh.cells if cb.type == "line"]
    )
    phys = (getattr(mesh, "cell_data_dict", {}) or {}).get("gmsh:physical", {})
    lineTag = (
        np.asarray(phys["line"]).ravel()
        if isinstance(phys, dict) and "line" in phys
        else None
    )
    if lineTag is not None and np.any(lineTag == tag):
        wallElems = lineElems[lineTag == tag]
    else:
        wallElems = lineElems
    pecNodes = np.unique(wallElems.ravel())
    X = nodeXY[tri]
    e = np.concatenate([X[:, 1] - X[:, 0], X[:, 2] - X[:, 1], X[:, 0] - X[:, 2]])
    minEdge = float(np.linalg.norm(e, axis=1).min())
    return nodeXY, tri, pecNodes, nodeXY[wallElems], minEdge


def assemble(nodeXY, tri):
    # Linear triangle: N_i = (a_i + b_i x + c_i y) / (2A), grad N_i = (b_i, c_i)/(2A)
    # with b = [y2-y3, y3-y1, y1-y2], c = [x3-x2, x1-x3, x2-x1],
    # A = |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)| / 2.
    # Stiffness: Ke_ij = int(grad N_i . grad N_j) = (b_i b_j + c_i c_j)/(4A).
    # Mass: Me_ij = int(N_i N_j) = A/12 * (1 + delta_ij),
    # i.e. A/12 * [[2,1,1],[1,2,1],[1,1,2]]. Assembled below in numpy
    X = nodeXY[tri]
    x1, y1 = X[:, 0, 0], X[:, 0, 1]
    x2, y2 = X[:, 1, 0], X[:, 1, 1]
    x3, y3 = X[:, 2, 0], X[:, 2, 1]
    b = np.stack([y2 - y3, y3 - y1, y1 - y2], axis=1)
    c = np.stack([x3 - x2, x1 - x3, x2 - x1], axis=1)
    area = 0.5 * np.abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    areaSafe = np.where(area == 0, 1.0, area)
    Ke = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (
        4.0 * areaSafe[:, None, None]
    )
    Me = area[:, None, None] * ((np.ones((3, 3)) + np.eye(3)) / 12.0)[None]
    K = np.zeros((len(nodeXY), len(nodeXY)))
    M = np.zeros((len(nodeXY), len(nodeXY)))
    I = np.broadcast_to(tri[:, :, None], (len(tri), 3, 3))
    J = np.broadcast_to(tri[:, None, :], (len(tri), 3, 3))
    np.add.at(K, (I, J), Ke)
    np.add.at(M, (I, J), Me)
    return K, M, b, c, area, X.mean(axis=1)


def solveMode(K, M, pecNodes, modeType, modeNum):
    # Generalized eigenproblem K v = kc^2 M v for -laplacian_t psi = kc^2 psi.
    # TM (psi = Ez): Dirichlet Ez = 0 on PEC, enforced by deleting PEC rows/cols.
    # TE (psi = Hz): Neumann dHz/dn = 0 is natural (zero boundary integral),
    # so the full matrices are used and the trivial constant mode (kc^2 ~ 0)
    # is skipped. Eigenvector normalized to max|psi| = 1.
    from scipy.linalg import eigh

    if modeType == "TM":
        free = np.setdiff1d(np.arange(K.shape[0]), np.unique(pecNodes))
        vals, vecs = eigh(
            K[np.ix_(free, free)],
            M[np.ix_(free, free)],
            subset_by_index=[0, modeNum - 1],
        )
        o = np.argsort(np.real(vals))
        psi = np.zeros(K.shape[0])
        psi[free] = np.real(vecs[:, o[modeNum - 1]])
        kc2 = float(np.real(vals)[o[modeNum - 1]])
    else:
        vals, vecs = eigh(K, M, subset_by_index=[0, modeNum])
        o = np.argsort(np.real(vals))
        vals, vecs = np.real(vals)[o], np.real(vecs)[:, o]
        keep = np.where(vals > 1e-3)[0][:modeNum]
        kc2 = float(vals[keep[-1]])
        psi = vecs[:, keep[-1]]
    m = np.max(np.abs(psi))
    return kc2, psi / m if m > 0 else psi


meshFile = sys.argv[1]
mode = sys.argv[2]
freq = float(sys.argv[3]) if len(sys.argv) > 3 else None
outDir = "results"

modeType = mode[:2].upper()
modeNum = int(mode[2:])

nodeXY, tri, pecNodes, pecWall, minEdge = loadMesh(meshFile)
K, M, b, c, area, triCent = assemble(nodeXY, tri)

# Free Space
epr = 1.0
mur = 1.0
eps = eps0 * epr
mu = mu0 * mur

# This doesn't work as expected for some cases
"""
if freq is None:
    kcTM, _ = solveMode(K, M, pecNodes, "TM", 1)
    kcTE, _ = solveMode(K, M, pecNodes, "TE", 1)
    fcTM = np.sqrt(max(kcTM, 0.0)) / (2 * np.pi * np.sqrt(eps * mu))
    fcTE = np.sqrt(max(kcTE, 0.0)) / (2 * np.pi * np.sqrt(eps * mu))
    freq = 2.0 * max(fcTM, fcTE)
    print("Used Frequency calculated:", freq / 1e6, "MHz")
"""
# Use 1Ghz directly as classroom comment
if freq is None:
    freq = 1e9

kc2, psi = solveMode(K, M, pecNodes, modeType, modeNum)
kc = np.sqrt(max(kc2, 0.0))
k = 2 * np.pi * freq * np.sqrt(eps * mu)
fc = kc / (2 * np.pi * np.sqrt(eps * mu))
s = k * k - kc * kc
beta = np.sqrt(s) if s >= 0 else 1j * np.sqrt(-s)
prop = s >= 0
omega = 2 * np.pi * freq

psiT = psi[tri]
gx = np.einsum("ei,ei->e", psiT, b) / (2 * np.where(area == 0, 1, area))
gy = np.einsum("ei,ei->e", psiT, c) / (2 * np.where(area == 0, 1, area))
bv = float(np.real(beta)) if prop else abs(complex(beta).imag)
if modeType == "TM":
    Ex = -(bv / kc**2) * gx
    Ey = -(bv / kc**2) * gy
    Hx = (omega * eps / kc**2) * gy
    Hy = -(omega * eps / kc**2) * gx
else:
    Ex = -(omega * mu / kc**2) * gy
    Ey = (omega * mu / kc**2) * gx
    Hx = -(bv / kc**2) * gx
    Hy = -(bv / kc**2) * gy
maxE = float(np.hypot(Ex, Ey).max(initial=0))
maxH = float(np.hypot(Hx, Hy).max(initial=0))
sclE = 0.9 * minEdge / maxE if maxE > 0 else 0.0
sclH = 0.9 * minEdge / maxH if maxH > 0 else 0.0
Ex, Ey = Ex * sclE, Ey * sclE
Hx, Hy = Hx * sclH, Hy * sclH
magE = np.hypot(Ex, Ey) / (0.9 * minEdge) if maxE > 0 else np.zeros_like(Ex)
magH = np.hypot(Hx, Hy) / (0.9 * minEdge) if maxH > 0 else np.zeros_like(Hx)

os.makedirs(outDir, exist_ok=True)
np.savetxt(
    outDir + "/longitudinal.dat",
    np.c_[nodeXY, psi],
    header="x y Re_psi(z=0)",
    comments="# ",
)
np.savetxt(
    outDir + "/transverse.dat",
    np.c_[triCent, Ex, Ey, Hx, Hy, magE, magH],
    header="xc yc Ex Ey Hx Hy magE magH",
    comments="# ",
)
f = open(outDir + "/boundary.dat", "w")
f.write("# PEC wall segments: x1 y1 x2 y2\n")
for (x1, y1), (x2, y2) in pecWall:
    f.write(f"{x1:.9e} {y1:.9e} {x2:.9e} {y2:.9e}\n")
f.close()

os.system("gnuplot plot.gp")

print("Calculated kc is :", kc, "rad/m")
print("Calculated fc is :", fc, "Hz")
print("Calculated beta is :", beta)
