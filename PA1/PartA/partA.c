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
// #define LivePlot
#define SaveGifPlot
// Also saves the figures in .pdf format
// #define TimeStampPlot
// Sim Constants
const double dt = 5E-13;
const int points = 2001;
const int decimationFactor = 30 * 4;
const int endTimeStep = 40000;

// Function which will handle all the sources
void excitationUpdate(Grid *g, int i, double dt) {
  g->Egrid[points / 2] =
      exp(-pow((i * dt - 5E-9) / 2E-9, 2)) * sin(2 * M_PI * 1E9 * i * dt);
}

int main() {
  Grid g;
  gplot g1;
  int decimationCounter = 0;

  // The domain was given to be from -3 to 3.
  allocateGrid(points, &g, 6.0);
  for (int i = 0; i < g.N; i++) {
    // Free Space
    g.epr[i] = 1.0;
    g.mur[i] = 1.0f;
  }

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
      "set yrange [-1.1:1.1];set y2range [-4E-3:4E-3]; set samples 800;"
      "set xrange [0:2001];");

  for (int i = 0; i < endTimeStep; i++) {
    forwardStep(&g, dt, 0);
    excitationUpdate(&g, i, dt);

#if defined(LivePlot) || defined(SaveGifPlot)
    if (!(decimationCounter++ % decimationFactor))
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric", "Magnetic", g.N);
#endif

#ifdef TimeStampPlot
    if (i == 13000) {
      gplotSetOption(&g1, "set output \"plot1.pdf\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric Field", "Magnetic Field",
                  g.N);
    }
    if (i == 37000) {
      gplotSetOption(&g1, "set output \"plot2.pdf\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric Field", "Magnetic Field",
                  g.N);
    }

#endif
    if (i == 22000) {
      printf("Calculated Velocity is :%lf m/s\n",
             returnVelocity(&g, 1746, dt, 2));

      printf("Calculated Impedance is :%lf Ohm\n",
             fabs(2 * g.Egrid[1620] / (g.Hgrid[1620] + g.Hgrid[1619])));
    }
  }
  gplotClose(&g1);
}
