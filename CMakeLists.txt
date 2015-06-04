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

cmake_minimum_required(VERSION 2.8)
project("efel")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)  

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")

exec_program(                                                                   
    "git"                                                                       
    ${CMAKE_CURRENT_SOURCE_DIR}                                                 
    ARGS "describe --tags --abbrev=0 --always"                                  
    OUTPUT_VARIABLE EFEL_MAINVERSION)                             
                                                                                
exec_program(                                                                   
    "git"                                                                       
    ${CMAKE_CURRENT_SOURCE_DIR}                                                 
    ARGS "describe --tags --abbrev=4 --dirty --always"                          
    OUTPUT_VARIABLE EFEL_VERSION)                                 
                                                                                
message("Working with efel version: ${EFEL_VERSION} (Main version: ${EFEL_MAINVERSION})")

if (BGQ)
    set(CMAKE_CXX_COMPILER mpixlcxx)
    set(CMAKE_C_COMPILER mpixlc)
    set(CMAKE_C_FLAGS "-qpic -qnostaticlink -O3 -qarch=qp -q64")
    set(CMAKE_CXX_FLAGS "-qpic -qnostaticlink -O3 -qarch=qp -q64")
elseif (LBGQ)
    set(CMAKE_CXX_COMPILER mpixlcxx)
    set(CMAKE_C_COMPILER mpixlc)
    set(CMAKE_C_FLAGS 
        "-qpic -qnohot -qstrict -qnostaticlink -O3 -qarch=qp -q64")
    set(CMAKE_CXX_FLAGS 
        "-qpic -qnohot -qstrict -qnostaticlink -O3 -qarch=qp -q64")
endif (BGQ)

add_subdirectory("efel/cppcore")
