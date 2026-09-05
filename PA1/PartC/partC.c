#include <stdio.h>
#define PLOT_IMPL
#include "../common/plotting.h"
#define SIM_IMPL
#include "../common/sim.h"
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

// Only one of these must be enabled
#define LivePlot
// #define SaveGifPlot
// Also saves the figures in .pdf format
// #define TimeStampPlot
// Sim Constants
const double dt = 5E-13;
const int points = 2001;
const int decimationFactor = 30 * 12;
const int endTimeStep = 20000;

// Function which will handle all the sources
void excitationUpdate(Grid *g, int i, double dt) {
  g->Egrid[points / 2] =
      (1 - exp(-dt * i / 1E-9)) * sin(2 * M_PI * 1E9 * i * dt);
}
void excitationUpdate2(Grid *g, int i, double dt, double phi) {
  g->Egrid[points / 2] =
      (1 - exp(-dt * i / 1E-9)) * sin((2 * M_PI * 1E9 * i * dt) + phi);
}
int main() {
  Grid g;
  Grid g2;
  gplot g1;
  int decimationCounter = 0;

  // Assuming the domain was given to be from -3 to 3.
  allocateGrid(points, &g, 6.0);
  for (int i = 0; i < g.N; i++) {
    // Free Space
    g.epr[i] = 1.0;
    g.mur[i] = 1.0f;
  }
  // Assuming the domain was given to be from -3 to 3.
  allocateGrid(points, &g2, 6.0);
  for (int i = 0; i < g2.N; i++) {
    // Free Space
    g2.epr[i] = 1.0;
    g2.mur[i] = 1.0f;
  }

  gplotInit(&g1);
  gplotSetOption(
#ifdef LivePlot
      &g1,
      "set term qt linewidth 1.5 size 1280,720;"
#endif
#ifdef SaveGifPlot
          & g1,
      "set term gif enhanced lw 2 size 720,720 animate delay 3 loop 0;"
#endif
#ifdef TimeStampPlot
          & g1,
      "set term pdfcairo lw 2 size 8in, 8in fontscale 0.8;"
#endif
      "set yrange [-1.1:1.1];set xrange [-1.1:1.1]; set samples 800; set "
      "xlabel \"E - along X\"; set ylabel \"E - along Y\"");

  double *Ey = (double *)calloc(endTimeStep, sizeof(double));
  double *Ex = (double *)calloc(endTimeStep, sizeof(double));
  char outputStringBuf[64];
  for (double phi = 0.0; phi < 2 * M_PI; phi += M_PI_4) {
    // Reset the Grid for every simulation
    decimationCounter = 0;
    memset(g.Egrid, 0, sizeof(double) * g.N);
    memset(g.Hgrid, 0, sizeof(double) * g.N);
    memset(g2.Egrid, 0, sizeof(double) * g.N);
    memset(g2.Hgrid, 0, sizeof(double) * g.N);
    memset(Ex, 0, sizeof(double) * endTimeStep);
    memset(Ey, 0, sizeof(double) * endTimeStep);

#if defined(SaveGifPlot)
    gplotSetOption(
        &g1,
        "set term gif enhanced lw 2 size 720,720 animate delay 3 loop 0;"
        "set yrange [-1.1:1.1];set xrange [-1.1:1.1]; set samples 800; set "
        "xlabel \"E - along X\"; set ylabel \"E - along Y\"");
    sprintf(outputStringBuf, "set output \"PlotPhi%.3lf.gif\"", phi);
    gplotSetOption(&g1, outputStringBuf);
#endif
#if defined(TimeStampPlot)
    sprintf(outputStringBuf, "set output \"PlotPhi%.3lf.pdf\"", phi);
    gplotSetOption(&g1, outputStringBuf);
#endif

    for (int i = 0; i < endTimeStep; i++) {
      forwardStep(&g, dt, 0);
      forwardStep(&g2, dt, 0);
      excitationUpdate(&g, i, dt);
      excitationUpdate2(&g2, i, dt, phi);
      Ey[i] = g.Egrid[1200];
      Ex[i] = g2.Egrid[1200];

#if defined(LivePlot) || defined(SaveGifPlot)
      if (!(decimationCounter++ % decimationFactor)) {
        // gplotPBDSEM(&g1, g2.Egrid, g2.Hgrid, "E", "H", g2.N);
        gplotXY(&g1, Ex, Ey, "Polarisation", i + 1);
        // printf("Plotting.\n");
      }
#endif
    }
#ifdef TimeStampPlot
    gplotXY(&g1, Ex, Ey, "polarisation", endTimeStep);
#endif
#ifdef SaveGifPlot
    gplotSetOption(&g1, "reset");
#endif
  }
  gplotClose(&g1);
}
