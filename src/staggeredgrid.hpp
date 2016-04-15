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
#ifndef STAGGEREDGRID_HPP
#define STAGGEREDGRID_HPP
#include <vector>
#include <functional>

#include "red3api.hpp"
#include "arrayops.hpp"
#include "array.hpp"

namespace red3
{

class RED3_API StaggeredGrid
{
public:
  StaggeredGrid(unsigned nnx, unsigned nny, unsigned nnz=1);
  virtual ~StaggeredGrid();

  double u(size_t i, size_t j, size_t k) { return m_u[INDEX(i,j,k,nu,nj,nk)]; }
  double v(size_t i, size_t j, size_t k) { return m_v[INDEX(j,i,k,nv,ni,nk)]; }
  double w(size_t i, size_t j, size_t k) { return m_w[INDEX(k,j,i,nw,nj,ni)]; }
  double p(size_t i, size_t j, size_t k) { return m_p[INDEX(i,j,k,nk,nj,nk)]; }
  void divg(Array<StaggeredGrid> &g);
  void dudx(Array<StaggeredGrid> &g);

  //void setU(std::function<double(double, double)> f);
  void setU(double(*f)(double x, double y));
  void setV(double(*f)(double x, double y));

  //void setU(std::function<double(double, double, double)> f);
  void setU(double(*f)(double x, double y, double z));
  void setV(double(*f)(double x, double y, double z));
  void setW(double(*f)(double x, double y, double z));

  double *allocateVariable();
  void divg(double *g);

  //double x(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double y(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double z(size_t i, size_t j, size_t k) { return m_u[i]; }

  //double xm(size_t i, size_t j, size_t k) { return m_xm[i]; }
  //double ym(size_t i, size_t j, size_t k) { return m_ym[i]; }
  //double zm(size_t i, size_t j, size_t k) { return m_zm[i]; }

  const unsigned ni;
  const unsigned nj;
  const unsigned nk;

  const unsigned nu;
  const unsigned nv;
  const unsigned nw;

  std::vector<double> x;
  std::vector<double> dx;
  std::vector<double> y;
  std::vector<double> dy;
  std::vector<double> z;
  std::vector<double> dz;
  std::vector<double> xm;
  std::vector<double> ym;
  std::vector<double> zm;

protected:
  unsigned m_nx;
  unsigned m_ny;
  unsigned m_nz;
  unsigned m_nu;
  unsigned m_nv;
  unsigned m_nw;
  unsigned m_ncells;

  double m_reynum;
  double m_dt;

  std::vector<double> m_x;
  std::vector<double> m_xm;
  std::vector<double> m_y;
  std::vector<double> m_ym;
  std::vector<double> m_z;
  std::vector<double> m_zm;

  double *m_u;
  double *m_v;
  double *m_w;
  double *m_p;

//private:
//  bool allocateSolution();
};

}

#endif