#pragma once
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  FILE *pipe;
} gplot;

void gplotInit(gplot *g);
void gplotSetOption(gplot *g, char *str);
void gplotPBD(gplot *g, double *data, char *title, int count);
void gplotPBDSEM(gplot *g, double *Edata, double *Hdata, char *title1,
                 char *title2, int count);
void gplotThree(gplot *g, double *d1, double *d2, double *d3, char *t1,
                char *t2, char *t3, int count);
void gplotXY(gplot *g, double *X, double *Y, char *title, int count);
void gplotClose(gplot *g);

#define PLOT_IMPL
#ifdef PLOT_IMPL

void gplotInit(gplot *g) { g->pipe = popen("gnuplot -persist\n", "w"); }

void gplotSetOption(gplot *g, char *str) { fprintf(g->pipe, "%s\n", str); }

void gplotPBD(gplot *g, double *data, char *title, int count) {
  fprintf(g->pipe,
          "plot \"-\" binary array=%d format='%%double' w l "
          "title \"%s\"\n",
          count, title);
  fwrite(data, sizeof(double), count, g->pipe);
  fflush(g->pipe);
}

void gplotPBDSEM(gplot *g, double *Edata, double *Hdata, char *title1,
                 char *title2, int count) {
  gplotSetOption(g, "set ytics nomirror; set y2tics");
  gplotSetOption(g, "set ylabel 'E'; set y2label 'H'");
  fprintf(g->pipe,
          "plot \"-\" binary array=%d format='%%double' w l "
          "smooth csplines title \"%s\" axes x1y1, "
          "\"-\" binary array=%d format='%%double' w l "
          "smooth csplines title \"%s\" axes x1y2\n",
          count, title1, count, title2);
  fwrite(Edata, sizeof(double), count, g->pipe);
  fwrite(Hdata, sizeof(double), count, g->pipe);
  fflush(g->pipe);
}

void gplotThree(gplot *g, double *d1, double *d2, double *d3, char *t1,
                char *t2, char *t3, int count) {
  fprintf(g->pipe,
          "plot \"-\" binary array=%d format='%%double' w l "
          "smooth csplines title \"%s\","
          "\"-\" binary array=%d format='%%double' w l "
          "smooth csplines title \"%s\","
          "\"-\" binary array=%d format='%%double' w l "
          "smooth csplines title \"%s\"\n",
          count, t1, count, t2, count, t3);
  fwrite(d1, sizeof(double), count, g->pipe);
  fwrite(d2, sizeof(double), count, g->pipe);
  fwrite(d3, sizeof(double), count, g->pipe);
}

void gplotXY(gplot *g, double *X, double *Y, char *title, int count) {
  double *buf = (double *)calloc(2 * count, sizeof(double));
  for (int j = 0; j < count; j++)
    buf[2 * j] = X[j];
  for (int j = 0; j < count; j++)
    buf[2 * j + 1] = Y[j];

  fprintf(g->pipe,
          "plot \"-\" binary record=%d format='%%double%%double' u 1:2 w l "
          "title \"%s\";\n",
          count, title);
  fwrite(buf, sizeof(double), 2 * count, g->pipe);
  free(buf);
}

void gplotClose(gplot *g) {
  fprintf(g->pipe, "exit\n");
  pclose(g->pipe);
  g->pipe = NULL;
}

#endif
