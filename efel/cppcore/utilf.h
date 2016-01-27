// Copyright © BBP/EPFL 2005-2011; All rights reserved. Do not distribute
// without further notice
/*----------------------------------------------------------------------------*/
/* SCCS Information: %W%    %G% */
/*----------------------------------------------------------------------------*/
#ifndef UTILF_H
#define UTILF_H

#include <math.h>

#ifndef PI
#define PI 3.14159265359
#endif

#define sq(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define POWER(x, y) exp((y) * log(x))

typedef double real;
typedef real *real1D;
typedef real **real2D;

#define realreadSI doublereadSI

extern real findmax(real2D u, int nx, int ny, int *imax, int *jmax);
extern real fmodmax(real2D u, int nx, int ny, int *imax, int *jmax);
extern real fmean(real2D u, real umean, int xs, int nx, int ny);
extern real findmin(real2D u, int nx, int ny, int *imin, int *jmin);
extern real sumfield(real2D u, int nx, int ny);
extern real findninf(real2D u, real2D v, int nx, int ny, int *imax, int *jmax);
extern real dp(real2D p, real2D cc, int nx, int ny);
/* compute the maximum relative error on divergence i.e.
   max(abs(div()/volume())) */
extern real rdivmax(real2D div, int nx, int ny, int *imax, int *jmax);
extern void copy(real2D out, real2D in, int nx, int ny);
extern void add(real2D a, real2D b, int nx, int ny);
extern void fill0(real2D u, int nx, int ny);

#endif
