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
// #define SaveGifPlot
// Also saves the figures in .pdf format
#define TimeStampPlot
// Sim Constants
const double dt = 5E-13;
// Choosen 6001 points because I need to divide the
// domain as 0 -- 2999 | 3000 | 3001 -- 6000
const int points = 6001;
const int decimationFactor = 30 * 4;
const int endTimeStep = 40000;

// Function which will handle all the sources
void excitationUpdate(Grid *g, int i, double dt) {
  g->Egrid[points / 2] +=
      exp(-pow((i * dt - 5E-9) / 2E-9, 2)) * sin(2 * M_PI * 1E9 * i * dt);
}

// Not handling s > e
double findMax(double *start, double *end) {
  double max = -INFINITY;
  while (start <= end) {
    max = (max < *start) ? *start : max;
    start++;
  }
  return max;
}

// Not handling s > e
double findMin(double *start, double *end) {
  double min = INFINITY;
  while (start <= end) {
    min = (min > *start) ? *start : min;
    start++;
  }
  return min;
}

int main() {
  Grid g;
  gplot g1;
  int decimationCounter = 0; // For plotting
  double incidentMax = 0, transmittedMax = 0, reflectedMax = 0;

  // The domain was given to be from -3 to 3.
  allocateGrid(points, &g, 6.0);
  for (int i = 0; i < g.N; i++) {
    // Free Space
    g.epr[i] = (i > 5000) ? 3.2 : 1.0;
    g.mur[i] = 1.0;
  }
  double *poyntingVec1 = calloc(endTimeStep, sizeof(double));
  double poyntingInt1 = 0, poyntingInt2 = 0;

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
      "set yrange [-5.1:5.1];set y2range [-2E-2:2E-2]; set samples 800;"
      "set xrange [0:6001]; set object 1 rect from graph 0.8333, graph 0 to "
      "graph 1, graph 1 behind fillcolor \"red\" fillstyle solid 0.3 noborder");

  for (int i = 0; i < endTimeStep; i++) {
    forwardStep(&g, dt, 0);
    excitationUpdate(&g, i, dt);

#if defined(LivePlot) || defined(SaveGifPlot)
    if (!(decimationCounter++ % decimationFactor))
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric", "Magnetic", g.N);
#endif

#ifdef TimeStampPlot
    if (i == 23000) {
      double tempE[g.N];
      double tempH[g.N];

      for (int i = 1; i < 2001; i++) {
        // E field is flipped
        tempE[3000 + i] = g.Egrid[3000 + i] - g.Egrid[3000 - i];
        tempH[3000 + i] = g.Hgrid[3000 + i] + g.Hgrid[3000 - i];
      }
      for (int i = 5001; i < 6001; i++) {
        tempE[i] = g.Egrid[i];
        tempH[i] = g.Hgrid[i];
      }
      // Reflected Wave
      gplotSetOption(&g1, "set output \"SeperateWave-TR.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Transmitted and Reflected Wave\"");
      gplotPBDSEM(&g1, tempE, tempH, "E", "H", g.N);
      for (int i = 3001; i < 5000; i++) {
        tempE[i] = g.Egrid[6000 - i];
        tempH[i] = -g.Hgrid[6000 - i];
      }
      // Incident wave
      gplotSetOption(&g1, "set output \"SeperateWave-TI.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Transmitted and Incident Wave\"");
      gplotPBDSEM(&g1, tempE, tempH, "E", "H", g.N);
      // Complete wave
      gplotSetOption(&g1, "set output \"SeperateWave-NW.pdf\"");
      gplotSetOption(&g1, "set xlabel \"Net Wave\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "E", "H", g.N);
    }

    if (i == 33333) {
      gplotSetOption(&g1, "set output \"plot1B.pdf\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric Field", "Magnetic Field",
                  g.N);
    }
    if (i == 39999) {
      gplotSetOption(&g1, "set output \"plot2B.pdf\"");
      gplotPBDSEM(&g1, g.Egrid, g.Hgrid, "Electric Field", "Magnetic Field",
                  g.N);
    }
#endif

    // Book-keeping for Part 4
    poyntingVec1[i] = -(g.Egrid[4996]) * (g.Hgrid[4996] + g.Hgrid[4997]) / 2;
    // printf("%lf\n", poyntingVec1[i]);
    poyntingInt1 -= (g.Egrid[4996]) * (g.Hgrid[4996] + g.Hgrid[4997]) / 2;
    poyntingInt2 -= (g.Egrid[5004]) * (g.Hgrid[5004] + g.Hgrid[5005]) / 2;

    if (i == 16000)
      incidentMax = (findMax(g.Egrid, g.Egrid + g.N - 1) >
                     fabs(findMin(g.Egrid, g.Egrid + g.N - 1)))
                        ? findMax(g.Egrid, g.Egrid + g.N - 1)
                        : findMin(g.Egrid, g.Egrid + g.N - 1);
    if (i == 33333) {
      printf("The velocity of the EM Wave in the Dielectric is: %lf\n",
             returnVelocity(&g, 5020, dt, 2));
      reflectedMax = (findMax(g.Egrid + 3000, g.Egrid + 4000) >
                      fabs(findMin(g.Egrid + 3000, g.Egrid + 4000)))
                         ? findMax(g.Egrid + 3000, g.Egrid + 4000)
                         : findMin(g.Egrid + 3000, g.Egrid + 4000);
    }
    if (i == 39999) {
      transmittedMax = (findMax(g.Egrid + 5001, g.Egrid + 6000) >
                        fabs(findMin(g.Egrid + 5001, g.Egrid + 6000)))
                           ? findMax(g.Egrid + 5001, g.Egrid + 6000)
                           : findMin(g.Egrid + 5001, g.Egrid + 6000);
    }
  }
  // printf("%lf %lf %lf\n", reflectedMax, transmittedMax, incidentMax);
  printf("Reflection Co-efficient is: %lf\n", (reflectedMax / incidentMax));
  printf("Transmission Co-efficient is: %lf\n", (transmittedMax / incidentMax));
  printf("Time Integral of Poynting Vector before and after the "
         "dielectric boundary: %lf %lf\n",
         poyntingInt1, poyntingInt2);
#ifdef TimeStampPlot
  gplotSetOption(&g1, "reset;set output \"PoyntingVec.pdf\"; set xlabel "
                      "\"TimeStep\"; set autoscale; set term pdfcairo lw 2 "
                      "size 8in, 6in fontscale 0.8");
  gplotPBD(&g1, poyntingVec1, "Poynting Vector before the dielectric",
           endTimeStep);
#endif
  gplotClose(&g1);
}
