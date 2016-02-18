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


#ifndef EFELLOGGER_H
#define EFELLOGGER_H

#include <iostream>
#include <fstream>
#include <vector>

class eFELLogger
{
public:
  eFELLogger(const string &outdir)
    : logging(false)
  {
    if (!outdir.empty()) {
      string filename = outdir + "/fllog.txt";
      logfile.open(filename.c_str(), std::fstream::out | std::fstream::app);
      logging = true;
    }
  }

  template<typename T>
  eFELLogger &operator<<(const std::vector<T> &v){
    if(logging){
      const size_t max = 10;
      for (size_t i = 0; i < v.size() && i < max; i++) {
        logfile << " " << v[i];
      }
      if (v.size() > max) {
        logfile << " ...";
      }
    }
    return *this;
  }

  template<typename T>
  eFELLogger &operator<<(const T value){
    if(logging){
      logfile << value;
    }
    return *this;
  }

  // deal with std::endl
  eFELLogger &operator <<(std::ostream& (*os)(std::ostream&)){
    if(logging){
      logfile << os;
    }
    return *this;
  }

private:
  bool logging;
  std::fstream logfile;
};

#endif
