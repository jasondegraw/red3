// Copyright (C) 2015-2017 Jason W. DeGraw
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
  StaggeredGrid(int ni, int nj, int nk=1, bool xperi=false);
  virtual ~StaggeredGrid()
  {}

  inline double u(int i, int j, int k) { return m_u[UINDEX(i, j, k, nu, nj, nk)]; }
  inline double v(int i, int j, int k) { return m_v[VINDEX(i, j, k, ni, nv, nk)]; }
  inline double w(int i, int j, int k) { return m_w[WINDEX(i, j, k, ni, nj, nw)]; }
  inline double p(int i, int j, int k) { return m_p[INDEX(i, j, k, nk, nj, nk)]; }
  //inline double var(int n, int i, int j, int k) { return m_var[n][INDEX(i, j, k, nk, nj, nk)]; }

  double northU(int i, int k) { return m_u_n[IKINDEX(i, k, nu, nk)]; }
  double southU(int i, int k) { return m_u_s[IKINDEX(i, k, nu, nk)]; }

  template <typename M> void divg(M &g)
  {
    int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    for(int k = 0; k < nk; k++) {
      for(int j = 0; j < nj; j++) {
        for(int i = 0; i < ni; i++) {
          g[i0] = (m_u[i1 + 1] - m_u[i1]) * m_rdx[i];
          i0++;
          i1++;
        }
        i1++;
      }
    }
    for(int k = 0; k < nk; k++) {
      for(int i = 0; i < ni; i++) {
        for(int j = 0; j < nj; j++) {
          g[j0] += (m_v[j1 + 1] - m_v[j1]) * m_rdy[j];
          j0++;
          j1++;
        }
        j1++;
      }
    }
    if(nw) {
      int k0 = 0, k1 = 0;
      for(int j = 0; j < nj; j++) {
        for(int i = 0; i < ni; i++) {
          for(int k = 0; k < nk; k++) {
            g[k0] += (m_w[k1 + 1] - m_w[k1]) * m_rdz[k];
            k0++;
            k1++;
          }
          k1++;
        }
      }
    }
  }
  void divg(Array<StaggeredGrid> &g);
  void divg(ArrayU<StaggeredGrid> &u, ArrayV<StaggeredGrid> &v, ArrayW<StaggeredGrid> &w, Array<StaggeredGrid> &g);
  void divg(ArrayU<StaggeredGrid> &u, ArrayV<StaggeredGrid> &v, Array<StaggeredGrid> &g);
  void dudx(Array<StaggeredGrid> &g);
  void dvdy(Array<StaggeredGrid> &g);

  //void setU(std::function<double(double, double)> f);
  void setU(double(*f)(double x, double y));
  void setV(double(*f)(double x, double y));

  //void setU(std::function<double(double, double, double)> f);
  void setU(double(*f)(double x, double y, double z));
  void setV(double(*f)(double x, double y, double z));
  void setW(double(*f)(double x, double y, double z));

  void setEastU(double(*f)(double y));
  void setWestU(double(*f)(double y));
  void setNorthU(double(*f)(double x));
  void setSouthU(double(*f)(double x));


  double *allocateVariable();
  void divg(double *g);

  //double x(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double y(size_t i, size_t j, size_t k) { return m_u[i]; }
  //double z(size_t i, size_t j, size_t k) { return m_u[i]; }

  //double xm(size_t i, size_t j, size_t k) { return m_xm[i]; }
  //double ym(size_t i, size_t j, size_t k) { return m_ym[i]; }
  //double zm(size_t i, size_t j, size_t k) { return m_zm[i]; }

  const int ni;
  const int nj;
  const int nk;

  const int nu;
  const int nv;
  const int nw;

  const bool xperi;

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
  int m_nx;
  int m_ny;
  int m_nz;
  int m_nu;
  int m_nv;
  int m_nw;
  int m_ncells;

  int m_dimension;

  std::vector<double> m_dx;
  std::vector<double> m_dy;
  std::vector<double> m_dz;

  std::vector<double> m_rdx;
  std::vector<double> m_rdy;
  std::vector<double> m_rdz;

  std::vector<double> m_x;
  std::vector<double> m_xm;
  std::vector<double> m_y;
  std::vector<double> m_ym;
  std::vector<double> m_z;
  std::vector<double> m_zm;

  std::unique_ptr<double[]> m_u;
  std::unique_ptr<double[]> m_v;
  std::unique_ptr<double[]> m_w;
  std::unique_ptr<double[]> m_p;

  double *m_u_n;
  double *m_u_s;

//private:
//  bool allocateSolution();
};

}

#endif