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

#ifndef __FILLFPTRTBL
#define __FILLFPTRTBL

#include "types.h"

#include "LibV1.h"
#include "LibV2.h"
#include "LibV3.h"
#include "LibV4.h"
#include "LibV5.h"

extern feature2function FptrTableV1;
extern feature2function FptrTableV2;
extern feature2function FptrTableV3;
extern feature2function FptrTableV4;
extern feature2function FptrTableV5;

int FillFptrTable();
#endif
