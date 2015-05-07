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


/* linear extrapolation in extra_velocity(), make_bc()
   and interface_tension() */
#define EXTRA_LINEAR

/* Upwind advection terms in momentum() */
#define UPWIND

/* Maximum number of multigrid V-cycles */
#define NCYCLEMAX 20


/* Maximum velocity in adaptative() */
#define MAXVEL 0.1

/* Multiplicating Factor for time step */
#define MULT 1.5

/* Dividing Factor for time step */
#define DIVIDE 2.0

/* Maximum number of V-cycle for adaptative() */
#define MAXCYCLE 5

/* Use axisymmetric smoothing of the markers */
#define AXISYMMETRIC_SMOOTHING

