#set page(paper: "a4", margin: (
  top: 1.5cm,
  left: 1.2cm,
  right: 1.2cm,
  bottom: 1cm,
))
#show underline: set underline(offset: 2pt)
#show link: set text(fill: blue)
#show link: underline
#set par(first-line-indent: 1.4em)
#align(center)[
  #show heading.where(level: 1): set text(size: 20pt)
  #show heading.where(level: 2): set text(size: 18pt)
  #set text(weight: "semibold", size: 14pt)
  = Waveguides and Antennas [EE3140]
  == Programming Assignment - 1
  #v(0.4em)
  #grid(
    column-gutter: 2em,
    columns: 2,
    [*EE24BTECH11201*], [*Shiva Varun Adicherla*],
  )
  #v(-0.4em)
  #line(length: 100%)
]

#set heading(numbering: "1.A")
#let pagebridge = align(right)[
  #v(1.2em)
  _(Continued in the next page)_
  #pagebreak()
]
#let to-sci(val, digits: 2) = {
  if val == 0 { return $0$ }
  let exponent = calc.floor(calc.log(calc.abs(val)))
  let mantissa = calc.round(val / calc.pow(10, exponent), digits: digits)
  $#mantissa times 10^(#exponent)$
}
#v(-1.1em)
// The 1st Question
= Part
My Assumptions will be Uniform waveguide along $z$ with PEC wall, homogeneous lossless dielectric ($sigma = 0$, real $mu, epsilon$, no sources). I'm going to use the same field
representation used in the lectures (Lecture 25). And I'm going to use *FEM* to solve this question. Since we'll be given a irregular mesh, and the equation which we are going
to solve came out to be of the form $nabla^2 psi = k_c^2 psi$, which is easy to solve.

== Transverse fields from longitudinal fields (Lecture 25)

The transverse parts of the two curl laws give
$
  bold(H)_perp = -1 slash (j omega mu)
  [nabla_perp times (E_z hat(bold(z))) + hat(bold(z)) times partial_z bold(E)_perp]
$
$
  bold(E)_perp = 1 slash (j omega epsilon)
  [nabla_perp times (H_z hat(bold(z))) + hat(bold(z)) times partial_z bold(H)_perp]
$
Eliminating $bold(H)_perp$ and using
$
  hat(bold(z)) times partial_z (hat(bold(z)) times partial_z bold(E)_perp)
  = -partial_z^2 bold(E)_perp
$
$ hat(bold(z)) times partial_z [nabla_perp times (E_z hat(bold(z)))]
= nabla_perp (partial_z E_z) $ yields, with $h^2 = omega^2 mu epsilon + gamma^2$,
$
  bold(E)_perp = -(j omega mu) slash h^2 nabla_perp times (H_z hat(bold(z)))
  - gamma slash h^2 nabla_perp E_z
$
$
  bold(H)_perp = (j omega epsilon) slash h^2 nabla_perp times (E_z hat(bold(z)))
  - gamma slash h^2 nabla_perp H_z
$
$
  nabla_perp times (psi hat(bold(z)))
  = (partial_y psi) hat(bold(x)) - (partial_x psi) hat(bold(y))
$
$ E_x = -(j omega mu) slash h^2 partial_y H_z - gamma slash h^2 partial_x E_z $
$ E_y = (j omega mu) slash h^2 partial_x H_z - gamma slash h^2 partial_y E_z $
$ H_x = (j omega epsilon) slash h^2 partial_y E_z - gamma slash h^2 partial_x H_z $
$ H_y = -(j omega epsilon) slash h^2 partial_x E_z - gamma slash h^2 partial_y H_z $
For $h != 0$, $E_z$ and $H_z$ cannot both vanish: TE is
$E_z = 0, H_z != 0$; TM is $H_z = 0, E_z != 0$. (TEM needs $h = 0$ and no
single hollow conductor supports it.) The above equations are taken from the lecture 25, and have been discussed thoroughly in the class.

== Helmholtz eigenproblem, cutoff, propagation constant

Following the convention of $exp(j omega t)$ for the forward moving wave, we have the following phasor equations.
$
  nabla times E & = -j omega mu H \
  nabla times H & = j omega epsilon E
$

Since the space inside the PEC boundaries is "free space" without any charge, we can write:
$
  nabla times nabla times E = -j omega mu (nabla times H) ==> nabla (0) - nabla^2 E = -j omega mu (j omega epsilon E) ==> nabla^2 E + omega^2 mu epsilon E = 0
$

Substituting only the longitudinal fields, we get: (A more rigorous derivation is in my rough note)

$ nabla_perp^2 E_z + h^2 E_z = 0, quad nabla_perp^2 H_z + h^2 H_z = 0 $
$-nabla_perp^2 psi = k_c^2 psi$ in $Omega$ with $psi = E_z$ (TM) or
$psi = H_z$ (TE) and $k_c^2 = h^2$. In all the general FEM Literature I've seen that the domain/mesh in which are solving the FEM is represented by $Omega$, and $k_c^2$ the eigen value,
and I'm doing the same here.

== PEC boundary conditions: derivation and enforcement

On a perfect conductor $n times E = 0$ and
$n dot B = 0$, with $n$ the outward cross-section normal the wall runs parallel to $z$, so $hat(z)$ is tangential to it. The $z$-component of $E$ is therefore tangential and must vanish:
$ E_z = 0 "on PEC (Dirichlet)". $
This is the TM condition. In the code it is enforced strongly, by all nodes on lines tagged `"pec"` are deleted from the system and $K_(F F) v = k_c^2 M_(F F) v$ is solved on the free nodes $F$, with
$psi = 0$ reinserted on the PEC nodes afterwards. Essentially removing the nodes from the solution/computation we decrease the degrees of freedom which is exactly this pinning of the field to 0 is doing.

For TE, $E_z = 0$ identically, so $n times E_perp = 0$ must constrain $H_z$ instead. From the TE reduction of the transverse formula,
$ E_perp = -(j omega mu) slash k_c^2 [(partial_y H_z) hat(x)) - (partial_x H_z) hat(y))] $
$ n_x partial_x H_z + n_y partial_y H_z = partial_n H_z = 0 #h(2em) "with" n = n_x hat(x) + n_y hat(y) $
$ partial_n H_z = 0 "on PEC (homogeneous Neumann)". $

Since the next step needs some FEM results, we'll go through them now. We have to solve the equation below:
$ nabla_perp^2 psi = - k_c^2 psi #h(1em) ("And as said before," psi "can be either" E "or" H "depending on mode") $
Now we multiply a basis function what has some properties which we can exploit down the line, like how we have lorentz's gauge condition helping us in simplifying EM equations in Antennas
(Magnetic Vector Potential A)
$
  integral_Omega w nabla_perp^2 psi dif Omega = -k_c^2 integral_Omega w psi dif Omega \
  "Using the integration by parts identity"\
  integral_Omega nabla w dot nabla psi dif Omega - integral_(partial Omega) w partial_n psi dif S = k_c^2 integral_Omega w psi dif Omega
$

Now if we come back and substitute $H$ instead of $psi$, and observe the equation, we are really only concerned for the middle term. But we chose the test/basis functions such that the
boundary integral vanishes i.e. $partial_n psi = 0$, so the code assembles the full $K, M$ ($K v = k_c^2 M v$) with no row elimination and simply discards the
trivial constant mode ($k_c approx 0$, non-physical for a hollow guide).

== Results of FEM simulation of a Rectangular waveguide (a = 2m, b = 1m)
My simulator just prints the data in the order of their cutoff frequencies since the *_eigh_* solver inside numpy directly returns an array in which eigen values are sorted, along with the corresponding
eigen vectors.

=== Results of TE Mode 3

#align(center, rect(
  ```
  ❯ python FEM.py rect.msh TE3

  Calculated kc is : 3.142447383217628 rad/m
  Calculated fc is : 149937011.10072228 Hz
  Calculated beta is : 20.721526489321995
  ```,
  outset: 0.2em,
  radius: 0.5em,
))
$
  h = k_c = sqrt(((m pi)/a)^2 + ((n pi)/b)^2) #h(2em) ("From the Lecture Notes")
$

Now, for TE3 in my simulator, we go through the modes one by one, $ "TE"_((1,0)) = pi/2 #h(1em) "TE"_((0,1)) = pi #h(1em) "TE"_((2,0)) = pi $

We now also notice that the cutoffs for TE2 and TE3 are going to be the same, and due to numerical precision errors accumilating differently for these different modes, we'll see either of them get
classified with higher or lower positions, and hence this ordering must not be expected always ($"TE"_(0,1) < "TE"_(2,0)$).

Coming back, we notice that our calculated $h = 3.14244738 "rad/"m$ which is close to what we are expecting. The error is:
$
  (delta h)/h = (3.14244738 - pi)/pi = 2.74 times 10^(-4)
$
#v(-1.5em)
#align(center, block(
  figure(
    image("PartA/results_TE3_rect/heatmap.pdf"),
    caption: [The Longitudinal Field (H) in this case, because TE mode],
  ),
  height: 1fr,
))
#align(center, block(
  figure(
    image("PartA/results_TE3_rect/transverse_H.pdf"),
    caption: [The Transverse H field in this Mode],
  ),
  height: 1fr,
))
#align(center, block(
  figure(
    image("PartA/results_TE3_rect/transverse_E.pdf"),
    caption: [The Transverse E field in this Mode],
  ),
  height: 1fr,
))

#v(1em)
#pagebreak()
=== Results of TM Mode 5
#align(center, rect(
  ```
  ❯ python FEM.py rect.msh TM5

  Calculated kc is : 7.03430513384043 rad/m
  Calculated fc is : 335630977.4894511 Hz
  Calculated beta is : 19.742724910406633
  ```,
  outset: 0.2em,
  radius: 0.5em,
))
$
  h = k_c = sqrt(((m pi)/a)^2 + ((n pi)/b)^2) #h(2em) ("From the Lecture Notes")
$

Now, for TM5 in my simulator, we go through the modes one by one,
$
  "TM"_((1,1)) = (pi sqrt(5))/2 #h(1em) "TM"_((2,1)) = pi sqrt(2) #h(1em) "TM"_((3,1)) = (pi sqrt(13))/2 #h(1em) "TM"_((1,2)) = (pi sqrt(17))/2 #h(1em) "TM"_((2,2)) = pi sqrt(5)
$

There's the next mode $"TM"_((4,1))$ that also has the same cutoff or eigen value, and again as mentioned the order is not concrete, and I've observed the output of the simulator and placed
$"TM"_((2,2))$ in the $5^("th")$ place.

Coming back, we notice that our calculated $h = 7.03430513 "rad/"m$ which is close to what we are expecting. The error is:
$
  (delta h)/h = (7.03430513 - pi sqrt(5))/(pi sqrt(5)) = 1.35 times 10^(-3)
$
#v(-1.5em)
#align(center, block(
  figure(
    image("PartA/results_TM5_rect/heatmap.pdf"),
    caption: [The Longitudinal Field (H) in this case, because TE mode],
  ),
  height: 30%,
))
#align(center, block(
  grid(
    columns: 2,
    figure(
      image("PartA/results_TM5_rect/transverse_H.pdf"),
      caption: [The Transverse H field in this Mode],
    ),
    figure(
      image("PartA/results_TM5_rect/transverse_E.pdf"),
      caption: [The Transverse E field in this Mode],
    ),
  ),
  height: 20%,
))

As we can see the Transverse $E$ and $H$ fields that I've plotted are always perpendicular to each other, like we expect in free space propagation.
The  results above were obtained for a pretty dense mesh, and I've observed that for lower density meshing, we had "clumping" of fields at certain nodes, which isn't expected since the Electric and
Magnetic field inside the waveguide were supposed to be smooth and continuous.

#pagebreak()
// The 2nd Question
= Part
In this question we are supposed to design a L-type matching networks for the following complex loads. Given $Z_0 = 75 Omega$ and $f_0 = 1 "Ghz"$.
== $Z_L = (150 + j 45) Omega$
To solve:
$
  z_l = (150+j 45)/75 & ==> 2 + j 3/5 = 2+ j 0.6#h(1em) "(Normalisation)"
$

#figure(
  image("PartB/150+j45_annotated.pdf", height: 51%),
  caption: [The work on ZY smith chart to design the L network],
)

As evident in the annotated smith chart, we are going to need a *reactive element* in parallel first, and then an additional element in series. The 2 possible solutions are the following, annotated
with green and blue sketches in the above chart. Since we are moving on Admittance Circles first, the *Parallel element comes first to the load*.
#v(-0.5em)
#columns(2, gutter: 0.3em)[
  On the Green line annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + green)[
    $
                     "Movement on Admittance Circle" \
      j omega L_("Parallel") = -1/(0.37 j) times Z_0 & ==> j 202 Omega \
              L_("Parallel") = 202/(2 pi times 10^9) & ==> 32.14 n H \
                      "Movement on Impedance Circle" \
         1/(j omega C_("Series")) = -1.1 j times Z_0 & ==> - j 82.5 Omega \
       C_("Series") = 1/(2 pi times 10^9 times 82.5) & ==> 1.92 p F
    $
  ])
  On the Cyan/Blue annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + blue)[
    $
                          "Movement on Admittance Circle" \
      1 / (j omega C_("Parallel")) = 1/(0.63 j) times Z_0 & ==> 119.04/j Omega \
        C_("Parallel") = 1/(2 pi times 10^9 times 119.04) & ==> 1.33 p F \
                           "Movement on Impedance Circle" \
                     j omega L_"Series" = 1.1 j times Z_0 & ==> 82.5j Omega \
                      L_"Series" = 82.5/(2 pi times 10^9) & ==> 13.13 n H
    $

  ])
]
#pagebreak()
#block(height: 1fr, width: 100%, figure(
  image("PartB/A-Blue.png"),
  caption: [The plot of the Blue annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a low pass, we can define the $3d b$ bandwidth to be $approx 1.9$Ghz
#block(height: 1fr, width: 100%, figure(
  image("PartB/A-Green.png"),
  caption: [The plot of the Green annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a high pass, we can't define the $3d b$ bandwidth, but $3d b$ point is $approx 0.4$Ghz
#v(3em)
#pagebreak()
== $Z_L = (30 - j 15) Omega$
To solve:
$
  z_l = (30-j 15)/75 & ==> 0.4 - j 0.2 #h(1em) "(Normalisation)"
$

#figure(
  image("PartB/30-j15_annotated.pdf", height: 57%),
  caption: [The work on ZY smith chart to design the L network],
)

As evident in the annotated smith chart, we are going to need a *reactive element* in parallel first, and then an additional element in series. The 2 possible solutions are the following, annotated
with green and orange sketches in the above chart. Here since we are travelling on the Impedance circles first, *Series element comes first to load*.
#v(-0.5em)
#columns(2, gutter: 0.3em)[
  On the Green line annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + green)[
    $
                    "Movement on Admittance Circle" \
      j omega L_("Parallel") = -1/(1.2 j) times Z_0 & ==> j 62.5 Omega \
            L_("Parallel") = 62.5/(2 pi times 10^9) & ==> 9.94 n H \
                     "Movement on Impedance Circle" \
        1/(j omega C_("Series")) = -0.3 j times Z_0 & ==> - j 22.5 Omega \
      C_("Series") = 1/(2 pi times 10^9 times 22.5) & ==> 7.07 p F
    $
  ])
  On the Orange annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + orange)[
    $
                         "Movement on Admittance Circle" \
      1 / (j omega C_("Parallel")) = 1/(1.2 j) times Z_0 & ==> 62.5/j Omega \
         C_("Parallel") = 1/(2 pi times 10^9 times 62.5) & ==> 2.54 p F \
                          "Movement on Impedance Circle" \
                    j omega L_"Series" = 0.7 j times Z_0 & ==> 52.5j Omega \
                     L_"Series" = 52.5/(2 pi times 10^9) & ==> 8.35 n H
    $

  ])
]
#pagebreak()
#block(height: 1fr, width: 100%, figure(
  image("PartB/B-Orange.png"),
  caption: [The plot of the Orange annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a low pass, we can define the $3d b$ bandwidth to be $approx 1$Ghz
#block(height: 1fr, width: 100%, figure(
  image("PartB/B-Green.png"),
  caption: [The plot of the Green annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a high pass, we can't define the $3d b$ bandwidth, but $3d b$ point is $approx 1$Ghz
#v(3em)
#pagebreak()

== $Z_L = (60+j 75) Omega$
To solve:
$
  z_l = (60+j 75)/75 & ==> 0.8 + j #h(1em) "(Normalisation)"
$

#figure(
  image("PartB/60+j75_annotated.pdf", height: 56%),
  caption: [The work on ZY smith chart to design the L network],
)

As evident in the annotated smith chart, we are going to need a *reactive element* in parallel first, and then an additional element in series. The 2 possible solutions are the following, annotated
with green and blue sketches in the above chart. The solution annotated with Blue travels on Impedance Circles first, so *Series element comes first to load for blue*, and similarly *Parallel element
comes first to load for green*.
#v(-0.5em)
#columns(2, gutter: 0.3em)[
  On the Blue line annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + blue)[
    $
                       "Movement on Admittance Circle" \
      1/(j omega C_("Parallel")) = 1/(0.5 j) times Z_0 & ==> 150/j Omega \
        C_("Parallel") = 1/(2 pi times 10^9 times 150) & ==> 1.06 p F \
                        "Movement on Impedance Circle" \
           1/(j omega C_("Series")) = -0.6 j times Z_0 & ==> - j 45 Omega \
           C_("Series") = 1/(2 pi times 10^9 times 45) & ==> 3.53 p F
    $
  ])
  On the Green annotated solution, we need:
  #v(-0.5em)
  #align(center, rect(outset: 0.3em, stroke: 1.8pt + green)[
    $
                         "Movement on Admittance Circle" \
      1 / (j omega C_("Parallel")) = 1/(0.1 j) times Z_0 & ==> 750/j Omega \
          C_("Parallel") = 1/(2 pi times 10^9 times 750) & ==> 0.212 p F \
                          "Movement on Impedance Circle" \
              1/(j omega C_"Series") = -1.02 j times Z_0 & ==> -76.5 j Omega \
             C_"Series" = 1/(2 pi times 10^9 times 76.5) & ==> 2.08 p F
    $

  ])
]
#pagebreak()
#block(height: 1fr, width: 100%, figure(
  image("PartB/C-Blue.png"),
  caption: [The plot of the Blue annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a high pass, we can't define the $3d b$ bandwidth, but $3d b$ point is $approx 0.1$Ghz
#block(height: 1fr, width: 100%, figure(
  image("PartB/C-Green.png"),
  caption: [The plot of the Green annotated Solution, Purple is the Gamma and Green is Transfer function (V)],
))
#v(1em)
As we can see, we have a good match at $f=1$Ghz as expected, and since the Transfer function looks like it's a high pass, we can't define the $3d b$ bandwidth, but $3d b$ point is $approx 0.3$Ghz
And the desmos graphs I've used to get these plots is available #link("https://www.desmos.com/calculator/n3dzb8yebg", "here").
#v(3em)
#pagebreak()


// The 3rd Question
= Part
==
As discussed in the class, if we assume that $omega$ is high enough to satisfy the conditions: $R << omega L$ and $G << omega C$, then we can consider that transmission line, a *_low loss_* transmission
line case, and the following quantities can be calculated by:

$
  Z_0 = sqrt((R+j omega L)/(G + j omega C)) approx sqrt(L/C) & ==> Z_0 = sqrt((250 times 10^(-9))/(100 times 10^(-12))) = 50 Omega \
  v_p = 1/sqrt(L C) & ==> v_p = 1/sqrt(25 times 10^(-18)) = 2 times 10^8 m s^(-1)
$

The above phase/propagation velocity is physical and valid, since it's lower than the intrinsic velocity of EM waves in Free Space / Vacuum ($c_0 approx 3 times 10^8 m s^(-1)$).

==
To determine the values of Reflection co-efficient $Gamma$ at the given 2 points, we:
$
              Gamma & = (Z_L - Z_0)/(Z_L + Z_0) \
    Gamma_"At Load" & = (100 - 50)/(100+50) ==> 1/3 \
    "In the Switch" & "is closed condition" \
  Gamma_"At Source" & = (25 - 50)/(25 + 50) ==> -1/3 \
    "In the Switch" & "is open condition" \
  Gamma_"At Source" & = (infinity - 50)/(infinity + 50) ==> 1 \
$

Some other important quantities that we are going to need are:
$
                   T_("Single Way") = d/v_p & => 20/(2 times 10^8) = 100 n s \
  V_("Initial") = 1 times Z_0 / (Z_0 + R_s) & => 50 / 75 = 2/3 V \
$

Also according to the question, we have the switch closed at $t=0s$ (Inferred it's was open before this), and then opens again at $t=250n s$.

In the case where we were calculating the $Gamma$ at the load and source, we are essentially seeing an impedance dis-continuity, and had to satisfy the general voltage and current equations
$V(z) = V_0^+ exp(-j beta z) + V_0^- exp(j beta z)$ and $I(z) = I_0^+ exp(-j beta z) - I_0^- exp(j beta z)$ continuity at the impedance mis-matched boundaries/z locations, along with satisfying $V(z) =
Z(z) I(z)$. When the switch opens, we'll have the current at the location of switch to be zero. To calculate the extra voltage pulse we get when opening the switch,
we solve the boundary condition in the following way.

$
  I(0, 250 n s) = 2/(3 Z_0) -2/(9 Z_0) -2/(27 Z_0) + V_2^+/Z_0 = 0 #h(3em) ("The 2/9 Voltage pulse is moving in opposite direction")
$

Then the Value of $V_2^+ = -10/27$

== Bounce Diagrams and Waveforms
#figure(
  image(
    "BoundsDiagram.png",
    width: 100%,
  ),
  caption: [The bounds diagram of the given scenario],
)

#import "@preview/simple-plot:1.1.0": line-plot, plot

#figure(
  plot(
    xmin: 0,
    xmax: 750,
    ymin: -0.2,
    ymax: 1.2,
    width: 17cm,
    height: 3.5cm,
    show-grid: true,
    xtick-step: 100,
    xlabel: [Time (ns)],
    ylabel: [Voltage (V)],

    line-plot(
      (
        (0, 0),
        (100, 0),
        (100, 0.899),
        (300, 0.899),
        (300, 0.790),
        (350, 0.790),
        (350, 0.296),
        (500, 0.296),
        (500, 0.263),
        (550, 0.263),
        (550, 0.098),
        (700, 0.098),
        (700, 0.087),
        (750, 0.087),
      ),
      stroke: blue + 1.2pt,
    ),
  ),
  caption: [The $V(t)$ observed at the Load],
)

#figure(
  plot(
    xmin: 0,
    xmax: 750,
    ymin: -0.2,
    ymax: 1.2,
    width: 17cm,
    height: 3.5cm,
    show-grid: true,
    xtick-step: 100,
    xlabel: [Time (ns)],
    ylabel: [Voltage (V)],

    line-plot(
      (
        (0, 0),
        (0, 0.6666),
        (200, 0.6666),
        (200, 0.814),
        (250, 0.814),
        (250, 0.444),
        (400, 0.444),
        (400, 0.395),
        (450, 0.395),
        (450, 0.148),
        (600, 0.148),
        (600, 0.131),
        (650, 0.131),
        (650, 0.049),
        (750, 0.049),
      ),
      stroke: blue + 1.2pt,
    ),
  ),
  caption: [The $V(t)$ observed at the Source],
)

The desmos versions of these graphs are available #link("https://www.desmos.com/calculator/nqwhf1hsvi", "at"), and we can also interactively view the voltages. However since we'll have
infinite reflections in this case, they've been truncated to about 5 in the calculations.

== FTDT Simulation
For this part, I've hacked/modified my previous FTDT simulator used in the previous assignment. All the changes have been documented in the source code as comments.
=== Lossless condition
The gist of this FTDT simulation is:
$
  (partial V)/(partial z) & = -L (partial I)/(partial t) - R I #h(3em) ("R is kept 0 in lossless condition") \
  (partial I)/(partial z) & = -C (partial V)/(partial t) - G V #h(3em) ("G is kept 0 in lossless condition") \
$
Which is solved by:
#pagebridge
#align(center)[
  #rect(
    ```c
    cCon = (dt / (grid->dx * Inductance));
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] -= cCon * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    cCon = (dt / (grid->dx * Capacitance));
    for (int j = 0; j < grid->N; j++) {
      grid->Egrid[j] -= cCon * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]

#figure(
  image(
    "PartC/Lossless_Voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

#let rows = csv("PartC/lossless.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.47,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}


=== Lossy condition
#v(0.5em)
#align(center)[
  #rect(
    ```c
    cCon = (dt / (2.0 * Inductance));
    double cCon1 = (1.0 - Resistance * cCon) / (1.0 + Resistance * cCon);
    double cCon2 = (dt / (Inductance * grid->dx)) / (1.0 + Resistance * cCon);
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] =
          cCon1 * grid->Hgrid[j] - cCon2 * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    cCon = (dt / (2.0 * Capacitance));
    cCon1 = (1.0 - Conductance * cCon) / (1.0 + Conductance * cCon);
    cCon2 = (dt / (Capacitance * grid->dx)) / (1.0 + Conductance * cCon);
    for (int j = 1; j < grid->N - 1; j++) {
      grid->Egrid[j] =
          cCon1 * grid->Egrid[j] - cCon2 * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]

#figure(
  image(
    "PartC/Lossy_voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

The Theoretical data in the following plots is still the lossless case, and is just there for Comparison. And as expected the peak/max voltages that we see in the simulation are lower than expected,
while keeping the exact same shape. This essentially means the transmission line is still dispersion-less.

#let rows = csv("PartC/lossy.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.47,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}


=== Dispersive - Dielectric Model
$
  epsilon_r (omega) = epsilon_oo + frac(epsilon_s - epsilon_oo, 1 + j omega tau) \
  C'(omega) = C_oo' + frac(Delta C', 1 + j omega tau), quad Delta C' = C_oo' (frac(epsilon_s - epsilon_oo, epsilon_oo)) \
  tau frac(partial P', partial t) + P'(t) = Delta C' V(t) \
  frac(partial V, partial z) = - L' frac(partial I, partial t) - R' I \
  frac(partial I, partial z) = - C_oo' frac(partial V, partial t) - frac(partial P', partial t) - G' V \
  v_(p, oo) = frac(1, sqrt(L' C_oo')), quad v_(p, 0) = frac(1, sqrt(L' (C_oo' + Delta C')))
$
Which is solved by:
#align(center)[
  #rect(
    ```c
    const double dC_p = Capacitance * (eps_s - eps_inf) / eps_inf;
    const double alpha = (2.0 * tau - dt) / (2.0 * tau + dt);
    const double beta_p = (2.0 * dt * dC_p) / (2.0 * tau + dt);
    const double cB = beta_p / (2.0 * dt);
    const double cP = (1.0 - alpha) / dt;

    cCon = (dt / (2.0 * Inductance));
    double cCon1 = (1.0 - Resistance * cCon) / (1.0 + Resistance * cCon);
    double cCon2 = (dt / (Inductance * grid->dx)) / (1.0 + Resistance * cCon);
    for (int j = 0; j < grid->N - 1; j++) {
      grid->Hgrid[j] =
          cCon1 * grid->Hgrid[j] - cCon2 * (grid->Egrid[j + 1] - grid->Egrid[j]);
    }

    memcpy(grid->Hgrid_prev, grid->Egrid_prev, sizeof(double) * grid->N);

    for (int j = 1; j < grid->N - 1; j++) {
      double dI = (grid->Hgrid[j] - grid->Hgrid[j - 1]) / grid->dx;
      double V_old = grid->Egrid[j];
      double P_old = grid->Hgrid_prev[j];

      double V_new = (V_old * (Capacitance / dt - cB - Conductance / 2.0) - dI +
                      cP * P_old) /
                      (Capacitance / dt + cB + Conductance / 2.0);

      grid->Egrid_prev[j] = alpha * P_old + beta_p * (V_new + V_old) / 2.0;
      grid->Egrid[j] = V_new;
    }
    ```,
    stroke: 0.6pt + black,
    outset: 0.6em,
    radius: 0.2em,
  ),
]
#v(-1em)
Where I'm using grid array structures grid->Egrid_prev and grid->Hgrid_prev as P(Polarisation) and P_old arrays

#figure(
  image(
    "PartC/Dielectric_Voltage.pdf",
    width: 65%,
  ),
  caption: [The Voltage across the transmission line at $65 n s$, the x-axis being the domain mesh points (0-20m)],
)

#let rows = csv("PartC/dielectric.csv")
#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(1))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (0, 0.6666),
          (200, 0.6666),
          (200, 0.814),
          (250, 0.814),
          (250, 0.444),
          (400, 0.444),
          (400, 0.395),
          (450, 0.395),
          (450, 0.148),
          (600, 0.148),
          (600, 0.131),
          (650, 0.131),
          (650, 0.049),
          (750, 0.049),
        ),
        stroke: blue + 1.2pt,
        label: [#text("Theoretical", fill: blue)],
        label-pos: 0.4,
      ),
      line-plot(
        points,
        stroke: red + 1.2pt,
        mark: none,
        label: [#text("Simulated", fill: red)],
      ),
    ),
    caption: [The $V(t)$ observed at the Source, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(2))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        (
          (0, 0),
          (100, 0),
          (100, 0.899),
          (300, 0.899),
          (300, 0.790),
          (350, 0.790),
          (350, 0.296),
          (500, 0.296),
          (500, 0.263),
          (550, 0.263),
          (550, 0.098),
          (700, 0.098),
          (700, 0.087),
          (750, 0.087),
        ),
        stroke: blue + 1.2pt,
        label: text("Theoretical", fill: blue),
        label-pos: 0.27,
      ),
      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Load, Comparison between Theoretical and Simulated],
  )
}

#{
  let points = rows
    .slice(1)
    .map(
      row => (float(row.at(0)) * 1E9, float(row.at(3))),
    )
  figure(
    plot(
      xmin: 0,
      xmax: 750,
      ymin: -0.2,
      ymax: 1.2,
      width: 17cm,
      height: 3.5cm,
      show-grid: true,
      xtick-step: 100,
      xlabel: [Time (ns)],
      ylabel: [Voltage (V)],

      line-plot(
        points,
        label: text("Simulated", fill: red),
        mark: none,
        stroke: red + 1.2pt,
      ),
    ),
    caption: [The $V(t)$ observed at the Mid point],
  )
}

And as we can see the output from the dielectric dispersion model is much more smeared, since the phase velocities for different frequencies aren't the same anymore, like in the case of Lossless or
Low-Loss case.

=== Comparison between the Theoretical/Analytical and Simulated results

As we can see from the above graphs, there's *severe ringing visible in the output Voltage and Current* which is the result of us *not having a smooth excitation source*, which will not contain the
infinite series of harmonics that can't exactly be represented in the finite discrete grid we have, and hence we can't expect the sharp rise and fall times in the voltage. #underline[A way to reduce
  the effect of this would be to use more finer meshing] but that is prohibitively expensive as we'd need much more compute time, and reality also won't let us have the ideal theoretical response we are
hoping for.

On the plus side, we see that there's a *very nice match between the expected amplitudes of Lossless Simulation and Analytical data*. And for the lossy case, #underline[since the expected loss was very low, we weren't able to observe much exponential decay, but the amplitude compared the lossless case was lower, while keeping the shape intact (So basically dispersion less), which was expected].

#v(1em)
#par(first-line-indent: 0em)[
  #text("GitHub Repository", weight: "bold", size: 16pt)
]
#v(0.2em)
The github repository can be accessed #link("https://github.com/ShivaVarunAdicherla/EE3140-PA", "here").
It will include a readme for this Programming assignment holding all the animated Gifs, and
the requested source code.

#par(first-line-indent: 0em)[
  #text("AI Disclosure", weight: "bold", size: 16pt)
]
#v(0.2em)
There was use of Generative AI (ChatGPT) to solve the Questions - 1 (Understanding the FEM and writing code for it, in python), and
Question - 3 (Idea to solve the Dielectric Model part), Code completely written by hand for the 3rd Question. All the use was related to coding and understanding the concept and problems, the report is
hand written.

#v(1em)
#line(length: 100%)
#align(center)[
  #text("END", weight: "bold", size: 14pt)
]

