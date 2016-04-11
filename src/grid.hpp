// Copyright (C) 2016 Jason W. DeGraw
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
#ifndef GRID_H
#define GRID_H

#include "red3api.hpp"
#include <vector>
#include <memory>

class RED3_API Grid1D
{
public:
  virtual ~Grid1D(){}
  virtual std::vector<double> grid() = 0;
  virtual std::vector<double> midgrid() = 0;
  virtual double delta(unsigned i) = 0;
  virtual unsigned n() = 0;
};

class RED3_API Uniform : public Grid1D
{
public:
  Uniform(unsigned n, double L=1, bool lengthIsDelta=false);
  virtual ~Uniform(){}
  virtual std::vector<double> grid();
  virtual std::vector<double> midgrid();
  virtual double delta(unsigned)
  {
    return m_dx;
  }
  virtual unsigned n()
  {
    return m_n;
  }

private:
  double m_L;
  unsigned m_n;
  double m_dx;
};

class RED3_API Grid2D
{
public:
  Grid2D(std::unique_ptr<Grid1D> x, std::unique_ptr<Grid1D> y);
};

#endif
