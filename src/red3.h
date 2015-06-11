// Copyright (C) 2015  Jason W. DeGraw
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef RED3_H
#define RED3_H

#include "red3api.h"

class RED3_API Grid1D
{
public:
  virtual bool setN(unsigned N) = 0;
  virtual bool setL(double L) = 0;
  virtual double x(unsigned i) = 0;
  virtual double xm(unsigned i) = 0;
};

class RED3_API Uniform
{
public:
  Uniform();
  virtual bool setN(unsigned n);
  virtual bool setL(double L);
  virtual double x(unsigned i);
  virtual double xm(unsigned i);
private:
  double m_L;
  unsigned m_n;
};

#endif
