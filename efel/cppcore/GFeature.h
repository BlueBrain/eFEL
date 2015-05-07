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

#ifndef __GFEATURE_H
#define __GFEATURE_H
int Initialize(const char *strDepFile, const char *outdir);
int setVersion(const char *strDepFile);
int setFeatureInt(const char *strName, int *A, int nValue);
int setFeatureDouble(const char *strName, double *A, int nValue);
int FeaturePrint(const char *strName);
int getFeatureInt(const char *strName, int **A);
int getFeatureDouble(const char *strName, double **A);
int setFeatureString(const char *key, const char *value);
int printFptr();
char *getgError();
double getDistance(const char *strName, double mean, double std);
#endif
