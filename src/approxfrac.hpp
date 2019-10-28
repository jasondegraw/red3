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
#ifndef RED3_APPROXFRAC_HPP
#define RED3_APPROXFRAC_HPP
#include "staggeredgrid.hpp"
#include "array.hpp"
#include <functional>
#include "red3api.hpp"
#include <Eigen/Sparse>

namespace red3 {
namespace approxfrac {

class RED3_API ViscousOperatorU : public StaggeredGrid::ArrayU
{
public:
  ViscousOperatorU(StaggeredGrid *parent, double dt, double Re);

  void compute(const StaggeredGrid::ArrayU &u);

  //StaggeredGrid::ArrayU last;
private:
  std::vector<double> m_ax;
  std::vector<double> m_bx;
  std::vector<double> m_cx;

  std::vector<double> m_ay;
  std::vector<double> m_by;
  std::vector<double> m_cy;

  std::vector<double> m_az;
  std::vector<double> m_bz;
  std::vector<double> m_cz;
};

class RED3_API ViscousOperatorV : public StaggeredGrid::ArrayV
{
public:
  ViscousOperatorV(StaggeredGrid *parent, double dt, double Re);

  void compute(const StaggeredGrid::ArrayV &u);

  //StaggeredGrid::ArrayU last;
private:
  std::vector<double> m_ax;
  std::vector<double> m_bx;
  std::vector<double> m_cx;

  std::vector<double> m_ay;
  std::vector<double> m_by;
  std::vector<double> m_cy;

  std::vector<double> m_az;
  std::vector<double> m_bz;
  std::vector<double> m_cz;
};

class RED3_API ViscousOperatorW : public StaggeredGrid::ArrayW
{
public:
  ViscousOperatorW(StaggeredGrid *parent, double dt, double Re);

  void compute(const StaggeredGrid::ArrayW &u);

  //StaggeredGrid::ArrayU last;
private:
  std::vector<double> m_ax;
  std::vector<double> m_bx;
  std::vector<double> m_cx;

  std::vector<double> m_ay;
  std::vector<double> m_by;
  std::vector<double> m_cy;

  std::vector<double> m_az;
  std::vector<double> m_bz;
  std::vector<double> m_cz;
};

class RED3_API IsothermalFlow : public StaggeredGrid
{
public:
  enum class BoundaryCondition {  };
  enum class Differencing {  };
  IsothermalFlow(double reynum, double dt, Grid1D& x, Grid1D& y, Grid1D& z=Grid1D::one(), bool xperi = false)
    : StaggeredGrid(x, y, z, xperi),
    g(this), //aw(this), an(this), as(this), af(this), ab(this), b(this), 
    reynum(reynum), dt(dt)
  {

  }


  void setupU();
  void setupV();
  void setupW();

  ArrayP g;
  const double reynum;
  const double dt;

  std::function<double(double)> A;

private:
  void convec(ArrayM &h1, ArrayM &h2, ArrayM &h3, index_t i1, index_t i2, index_t j1, index_t j2, index_t k1, index_t k2);

};

}
}

#endif