#include <stdio.h>
#define PLOT_IMPL
#include "../common/plotting.h"
#define SIM_IMPL
#include "../common/sim.h"
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

// I'm assuming that Egrid will now be V or Potential, and Hgrid will be I or
// Current

// Only one of these must be enabled
// #define LivePlot
// #define SaveGifPlot
// Also saves the figures in .pdf format
#define TimeStampPlot
// Sim Constants
const double dt = 5E-12;
const int points = 2001;
const int decimationFactor = 30 * 4;
const int endTimeStep = 150000;
const double Inductance = 250E-9;
const double Capacitance = 100E-12;
const double Resistance = 0.15;
const double Conductance = 1E-7;

void forwardStepTransmission(Grid *grid, double dt, int mode, int i) {
  double cCon;

#ifdef LOSSLESS
  cCon = (dt / (grid->dx * Inductance));
  for (int j = 0; j < grid->N - 1; j++) {
    grid->Hgrid[j] -= cCon * (grid->Egrid[j + 1] - grid->Egrid[j]);
  }

  cCon = (dt / (grid->dx * Capacitance));
  for (int j = 1; j < grid->N - 1; j++) {
    grid->Egrid[j] -= cCon * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
  }
#endif

#ifdef LOSSY
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
#endif
// AI's Idea, My initial idea was to decompose the step into different
// sinusoids, simulate it for all those sinusoids, then add them up. But AI's
// Idea was better. Debye ADE-FDTD: q' = C_p V + P with tau dP/dt + P = dC_p
// V, dC_p = C_p (eps_s - eps_inf) / eps_inf so the DC limit is C_p + dC_p.
#ifdef DIELEC
  const double eps_inf = 2.25;
  const double eps_s = 2.56;
  const double tau = 5E-9;
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
#endif

  // Boundary Conditions
  cCon = (Capacitance * grid->dx / (2.0 * dt));
  const double Gs_half = 1.0 / (2.0 * 25.0);  // 0.02 S
  const double GL_half = 1.0 / (2.0 * 100.0); // 0.005 S

  if (mode == 0) {
    if (i * dt < 250E-9) {
      grid->Egrid[0] =
          (1.0 / 25.0 + (cCon - Gs_half) * grid->Egrid[0] - grid->Hgrid[0]) /
          (cCon + Gs_half);
    } else {
      grid->Egrid[0] -= grid->Hgrid[0] / cCon;
    }

    grid->Egrid[grid->N - 1] = (grid->Hgrid[grid->N - 2] +
                                (cCon - GL_half) * grid->Egrid[grid->N - 1]) /
                               (cCon + GL_half);
  }
}

int main() {
  Grid g;
  gplot g1;
  int decimationCounter = 0;

  allocateGrid(points, &g, 20.0);
  for (int i = 0; i < g.N; i++) {
    // Free Space
    g.epr[i] = 1.0;
    g.mur[i] = 1.0f;
  }

  FILE *simData = fopen("./data.csv", "w");

  gplotInit(&g1);
  gplotSetOption(
#ifdef LivePlot
      &g1,
      "set term qt linewidth 1.5 size 1280,720;"
#endif
#ifdef SaveGifPlot
          & g1,
      "set term gif enhanced lw 2 size 720,480 animate delay 3 loop 0; set "
      "output \"plot.gif\";"
#endif
#ifdef TimeStampPlot
          & g1,
      "set term pdfcairo lw 2 size 8in, 4in fontscale 0.8;"
#endif
      "set yrange [-0.1:1.2]; set y2range [-0.1: 0.1]; set samples 800;"
      "set xrange [0:2001];");

  for (int i = 0; i < endTimeStep; i++) {
    forwardStepTransmission(&g, dt, 0, i);
    // excitationUpdate(&g, i, dt);

#if defined(LivePlot) || defined(SaveGifPlot)
    if (!(decimationCounter++ % decimationFactor))
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Voltage", "Current", g.N);
#endif

#ifdef TimeStampPlot
    if (i == 13000) {
      gplotSetOption(&g1, "set output \"plot1.pdf\"");
      gplotPBD(&g1, g.Egrid, "Voltage", g.N);
    }
#endif
    if (!(i % 1000))
      fprintf(simData, "%g,%g,%g,%g\n", i * dt, g.Egrid[1], g.Egrid[g.N - 1],
              g.Egrid[(g.N - 1) / 2]);
  }
  fclose(simData);
  gplotClose(&g1);
}
