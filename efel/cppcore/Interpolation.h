/* Copyright (c) 2015, EPFL/Blue Brain Project                                   
 *                                                                               
 * This file is part of eFEL <https://github.com/BlueBrain/eFEL>                 
 *                                                                               
 * This library is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License version 3.0 as published   
 * by the Free Software Foundation.                                              
 *                                                                               
 * This library is distributed in the hope that it will be useful, but WITHOUT   
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more 
 * details.                                                                      
 *                                                                               
 * You should have received a copy of the GNU Lesser General Public License      
 * along with this library; if not, write to the Free Software Foundation, Inc., 
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                   
 */        

#ifndef __INTERPOLATION_H                                                                  
#define __INTERPOLATION_H                                                                  
#define POLYMAX 100
#define SP_NATURAL 1
#define SP_DERIVATIVE 2
#define SP_PERIODIC 3

#define splint(poly, x) \
  (poly.a + (x) * (poly.b + (x) * (poly.c + (x) * poly.d)))
#define splint1(poly, x) (poly.b + (x) * (2. * poly.c + (x) * 3. * poly.d))
#define splint2(poly, x) (2. * poly.c + (x) * 6. * poly.d)

typedef struct {
  real a, b, c, d;
} polynom3;

extern real polynomial(real *xp, real *yp, real *z, real x, real y);

extern void CubicPolynom(real x0, real y0, real x1, real y1, real x2, real y2,
                         real x3, real y3, real *ex, real *fx, real *gx,
                         real *ey, real *fy, real *gy, real *rd2);

extern void spline(real1D x, real1D y, int n, real yp1, real ypn, real1D d2y,
                   polynom3 *poly, int flag);

extern void tridag(real1D a, real1D b, real1D c, real1D r, real1D u, int n);

extern void Pspline(real1D x, real1D y, real1D d2y, real lp, int n,
                    polynom3 *poly);

extern real splint_find(real1D xa, polynom3 *poly, int n, real x);

extern int distmin(polynom3 px, polynom3 py, real xp, real yp, real *tmin,
                   real t1, real t2, real prec);

extern int polyroot(polynom3 poly, real x1, real x2, real *x, real prec);

extern int polyroots(polynom3 poly, real x1, real x2, real *x, real y,
                     real prec);

extern void bcucof(real2D f, int i, int j, real *c, int nx, int ny);

extern real bcuint(real t, real u, real *c);

extern real polydistmin(polynom3 px1, polynom3 px2, polynom3 py1, polynom3 py2,
                        real *t1, real *t2, real ftol);

extern real trapzd(real (*func)(polynom3, polynom3, real), polynom3 px,
                   polynom3 py, real a, real b, int n);

extern real qtrap(real (*func)(polynom3, polynom3, real), polynom3 px,
                  polynom3 py, real a, real b, real eps);

#endif
