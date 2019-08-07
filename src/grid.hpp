// Copyright (C) 2015-2019 Jason W. DeGraw
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
#ifndef RED3_GRID_HPP
#define RED3_GRID_HPP

#include "red3.hpp"
#include <vector>
#include <memory>

namespace red3 {

class RED3_API Grid1D
{
public:
  enum class Generator {UniformN, UniformDelta};

  Grid1D(index_t N);

  virtual ~Grid1D(){}

  std::vector<double> grid() const
  {
    return m_x;
  }

  std::vector<double> midgrid() const
  {
    return m_xm;
  }

  std::vector<double> deltas() const
  {
    return m_dx;
  }

  double &operator[](index_t i)
  {
    return m_x[i];
  }

  double &grid(index_t i)
  {
    return m_x[i];
  }

  double &midgrid(index_t i)
  {
    return m_x[i];
  }

  double &m(index_t i)
  {
    return m_x[i];
  }

  double &delta(index_t i)
  {
    return m_dx[i];
  }

  double &delta0()
  {
    return m_dx[0];
  }

  double &deltaN()
  {
    return m_dx.back();
  }

  double operator[](index_t i) const
  {
    return m_x[i];
  }

  double grid(index_t i) const
  {
    return m_x[i];
  }

  double midgrid(index_t i) const
  {
    return m_x[i];
  }

  double m(index_t i) const
  {
    return m_x[i];
  }

  double delta(index_t i) const
  {
    return m_dx[i];
  }

  double delta0() const
  {
    return m_dx[0];
  }

  double deltaN() const
  {
    return m_dx.back();
  }

  index_t size() const
  {
    return N;
  }

  const Generator generator;
  const index_t N;
  const bool uniform;

private:
  std::vector<double> m_x;
  std::vector<double> m_xm;
  std::vector<double> m_dx;
};

class RED3_API Grid2D
{
public:
  Grid2D(std::unique_ptr<Grid1D> x, std::unique_ptr<Grid1D> y);
};

}

#endif
