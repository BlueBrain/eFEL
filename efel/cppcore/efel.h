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

#ifndef __EFEL_H                                                       
#define __EFEL_H 
#define FEATURELIB_API
extern "C" {
FEATURELIB_API int Initialize(const char *strDepFile, const char *outdir);
FEATURELIB_API int setVersion(const char *strDepFile);
FEATURELIB_API int setFeatureInt(const char *strName, int *A, unsigned nValue);
FEATURELIB_API int setFeatureDouble(const char *strName, double *A, unsigned nValue);
FEATURELIB_API int getTotalIntData();
FEATURELIB_API int getTotalDoubleData();
FEATURELIB_API int FeaturePrint(const char *strName);
FEATURELIB_API int getFeatureInt(const char *strName, int **A);
FEATURELIB_API int getFeatureDouble(const char *strName, double **A);
FEATURELIB_API int setFeatureString(const char *key, const char *value);
FEATURELIB_API int getFeatureString(const char *key, char **value);
FEATURELIB_API int printFptr();
FEATURELIB_API char *getgError();
FEATURELIB_API double getDistance(const char *strName, double mean, double std, bool trace_check);
}
#endif
