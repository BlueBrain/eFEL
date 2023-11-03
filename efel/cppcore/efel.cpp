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

#include "efel.h"
#include "cfeature.h"

cFeature *pFeature = NULL;

int Initialize(const char *strDepFile, const char *outdir) {
  if (pFeature != NULL) {
    delete pFeature;
  }

  pFeature = new cFeature(string(strDepFile), string(outdir));
  if (pFeature == NULL) {
    return -1;
  } else {
    return 1;
  }
}

char *getgError() {
  string error = GErrorStr + pFeature->getGError();
  GErrorStr.clear();
  return (char *)error.c_str();
}
