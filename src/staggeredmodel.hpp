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
#include <vector>

#include "red3api.hpp"
#include "arrayops.hpp"

namespace red3
{

class RED3_API StaggeredGrid
{
public:
  StaggeredGrid(unsigned nnx, unsigned nny, unsigned nnz=1);
  virtual ~StaggeredGrid();

  double u(size_t i, size_t j, size_t k) { return m_u[i]; }
  double v(size_t i, size_t j, size_t k) { return m_v[i]; }
  double w(size_t i, size_t j, size_t k) { return m_w[i]; }
  double p(size_t i, size_t j, size_t k) { return m_p[i]; }

  //double x(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double y(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double z(size_t i, size_t j, size_t k) { return m_u[i]; }

  //double xm(size_t i, size_t j, size_t k) { return m_xm[i]; }
  //double ym(size_t i, size_t j, size_t k) { return m_ym[i]; }
  //double zm(size_t i, size_t j, size_t k) { return m_zm[i]; }

  unsigned nx;
  unsigned ny;
  unsigned nz;

  unsigned nu;
  unsigned nv;
  unsigned nw;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
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
//  double *m_x;
//  double *m_xm;
//  double *m_y;
//  double *m_ym;
//  double *m_z;
//  double *m_zm;
  double *m_u;
  double *m_v;
  double *m_w;
  double *m_p;

private:
  bool allocateSolution();
};

}
