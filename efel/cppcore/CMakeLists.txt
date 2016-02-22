# Copyright (c) 2015, EPFL/Blue Brain Project                                   
#                                                                               
# This file is part of eFEL <https://github.com/BlueBrain/eFEL>                 
#                                                                                
# This library is free software; you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License version 3.0 as published   
# by the Free Software Foundation.                                              
#                                                                                
# This library is distributed in the hope that it will be useful, but WITHOUT   
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more 
# details.                                                                      
#                                                                                
# You should have received a copy of the GNU Lesser General Public License      
# along with this library; if not, write to the Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                   

cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE Debug)

set(FEATURESRCS Utils.cpp LibV1.cpp LibV2.cpp LibV3.cpp LibV4.cpp LibV5.cpp
    FillFptrTable.cpp DependencyTree.cpp efel.cpp cfeature.cpp
    mapoperations.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

add_library(efelStatic ${FEATURESRCS})
set_target_properties(efelStatic PROPERTIES OUTPUT_NAME efel)
install(TARGETS efelStatic ARCHIVE DESTINATION lib)

add_library(efel SHARED ${FEATURESRCS})
install(TARGETS efel LIBRARY DESTINATION lib)

install(FILES efel.h cfeature.h FillFptrTable.h LibV1.h LibV2.h LibV3.h
    LibV4.h LibV5.h mapoperations.h Utils.h DependencyTree.h eFELLogger.h
    types.h
    DESTINATION include)
