#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const double mu0 = 4 * M_PI * 1E-7;
const double ep0 = 8.854187817E-12;
const double c0 = 299792458.0;

typedef struct {
  double *Egrid;
  double *Egrid_prev;
  double *Hgrid;
  double *Hgrid_prev;
  int N;
  double dx;
  double *epr;
  double *mur;
} Grid;

void allocateGrid(int N, Grid *grid, double domain);
void forwardStep(Grid *grid, double dt, int mode);
double returnVelocity(Grid *g, int i, double dt, int window);

#define SIM_IMPL
#ifdef SIM_IMPL

void allocateGrid(int N, Grid *grid, double domain) {
  grid->Egrid = (double *)calloc(N, sizeof(double));
  grid->Egrid_prev = (double *)calloc(N, sizeof(double));
  grid->Hgrid = (double *)calloc(N, sizeof(double));
  grid->Hgrid_prev = (double *)calloc(N, sizeof(double));
  grid->N = N;
  grid->dx = (domain / (N - 1));
  grid->epr = (double *)calloc(N, sizeof(double));
  grid->mur = (double *)calloc(N, sizeof(double));
}

// Mode: 0 Both PEC, 1 Left is PEC, 2 Right is PEC, 3 None is PEC
void forwardStep(Grid *grid, double dt, int mode) {
  memcpy(grid->Egrid_prev, grid->Egrid, sizeof(double) * grid->N);
  memcpy(grid->Hgrid_prev, grid->Hgrid, sizeof(double) * grid->N);

  double cCon = (dt / (grid->dx * mu0));
  for (int j = 0; j < grid->N - 1; j++) {
    grid->Hgrid[j] +=
        (cCon / grid->mur[j]) * (grid->Egrid[j + 1] - grid->Egrid[j]);
  }

  cCon = (dt / (grid->dx * ep0));
  for (int j = 1; j < grid->N - 1; j++) {
    grid->Egrid[j] +=
        (cCon / grid->epr[j]) * (grid->Hgrid[j] - grid->Hgrid[j - 1]);
  }

  if (mode == 0) {
    grid->Egrid[0] = 0;
    grid->Egrid[grid->N - 1] = 0;
  }
}

double returnVelocity(Grid *g, int i, double dt, int window) {
  if (!g || i < 0 || i >= g->N || dt <= 0.0) {
    return 1.0;
  }

  double target = g->Egrid[i];
  int start = (i - window < 0) ? 0 : (i - window);
  int end = (i + window >= g->N - 1) ? (g->N - 2) : (i + window);

  double min_dist = -1.0;
  double best_prev_idx = 0.0;
  int found = 0;

  for (int j = start; j <= end; j++) {
    double e1 = g->Egrid_prev[j];
    double e2 = g->Egrid_prev[j + 1];

    if ((e1 <= target && target <= e2) || (e2 <= target && target <= e1)) {
      double alpha = 0.0;
      double diff = e2 - e1;
      if (fabs(diff) > 1e-13) {
        alpha = (target - e1) / diff;
      }

      double prev_idx = (double)j + alpha;
      double dist = fabs((double)i - prev_idx);

      if (!found || dist < min_dist) {
        min_dist = dist;
        best_prev_idx = prev_idx;
        found = 1;
      }
    }
  }

  if (!found) {
    return 1.0;
  }

  double dz = ((double)i - best_prev_idx) * g->dx;
  // printf("%e %e\n", dz, dt);
  return dz / dt;
}

#define EPS_ENERGY 1e-15 // Energy cutoff to prevent 0/0 division

#endif
