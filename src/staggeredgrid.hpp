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
#include <algorithm>
#include <memory>

#include "red3api.hpp"
#include "array.hpp"
#include "util.hpp"
#include "grid.hpp"

namespace red3
{

class RED3_API StaggeredGrid
{
public:
  using Array = ChildArray<StaggeredGrid>;
  using ArrayU = ChildArrayU<StaggeredGrid>;
  using ArrayV = ChildArrayV<StaggeredGrid>;
  using ArrayW = ChildArrayW<StaggeredGrid>;

  StaggeredGrid(int ni, int nj, int nk=1, bool xperi=false)
    : ni(ni), nj(nj), nk(std::max(nk,1)), nu(ni + 1), nv(nj + 1), nw(nk == 1 ? 1 : nk + 1), xperi(xperi)
  {
    // Check the inputs
    if (ni < 2) {
      fatal("ni must be greater than 1");
    }
    if (nj < 2) {
      fatal("nj must be greater than 1");
    }

    m_ncells = ni*nj*nk;

    int nuvw = ni + nj + 2;
    if (nk > 1) {
      m_nu = nu*nj*nk;
      m_nv = ni*nv*nk;
      m_nw = ni*nj*nw;
      w = wArray();
    }
    else {
      m_nu = nu*nj;
      m_nv = ni*nv;
      m_nw = 0;
    }
    u = uArray();
    v = vArray();
    p = array();
    m_u = std::make_unique<double[]>(m_nu); // (double*)callocate(m_nu, sizeof(double), "u velocity");
    m_u_n = (double*)callocate(ni*nk, sizeof(double), "north u velocity");
    m_u_s = (double*)callocate(ni*nk, sizeof(double), "south u velocity");
    m_v = std::make_unique<double[]>(m_nv); //(double*)callocate(m_nv, sizeof(double), "v velocity");
                                            //m_x = (double*)callocate(nu, sizeof(double), "x grid");
                                            //m_y = (double*)callocate(nv, sizeof(double), "y grid");
                                            //m_xm = (double*)callocate(nx, sizeof(double), "x mid-grid");
                                            //m_ym = (double*)callocate(ny, sizeof(double), "y mid-grid");
    m_p = std::make_unique<double[]>(m_ncells); //(double*)callocate(m_ncells, sizeof(double), "pressure");

    // Make a simple grid
    Uniform xg(nu);
    x = xg.grid();
    xm = xg.midgrid();
    dx = xg.deltas();
    double rdx = 1.0 / xg.delta0();
    m_rdx = std::vector<double>(xg.n(), rdx);
    Uniform yg(nv);
    y = yg.grid();
    ym = yg.midgrid();
    dy = yg.deltas();
    double rdy = 1.0 / yg.delta0();
    m_rdy = std::vector<double>(yg.n(), rdy);
    if (nk > 1) {
      Uniform zg(nv);
      z = zg.grid();
      zm = zg.midgrid();
      dz = zg.deltas();
      double rdz = 1.0 / zg.delta0();
      m_rdz = std::vector<double>(zg.n(), rdz);
    } else {
      z = { -0.5, 0.5 };
      zm = { 0.0 };
      dz = { 1.0 };
      m_rdz = { 1.0 };
    }
    
  }

  virtual ~StaggeredGrid()
  {}

  //inline double u(int i, int j, int k) { return m_u[UINDEX(i, j, k, nu, nj, nk)]; }
  //inline double v(int i, int j, int k) { return m_v[VINDEX(i, j, k, ni, nv, nk)]; }
  //inline double w(int i, int j, int k) { return m_w[WINDEX(i, j, k, ni, nj, nw)]; }
  //inline double p(int i, int j, int k) { return m_p[INDEX(i, j, k, nk, nj, nk)]; }
  //inline double var(int n, int i, int j, int k) { return m_var[n][INDEX(i, j, k, nk, nj, nk)]; }

  //inline double u(int ijk) { return m_u[ijk]; }
  //inline double v(int ijk) { return m_v[ijk]; }
  //inline double w(int ijk) { return m_w[ijk]; }
  //inline double p(int ijk) { return m_p[ijk]; }

  double northU(int i, int k) { return m_u_n[IKINDEX(i, k, nu, nk)]; }
  double southU(int i, int k) { return m_u_s[IKINDEX(i, k, nu, nk)]; }

  double eastU(int j, int k) { return m_u[UINDEX(0, j, k, nu, nj, nk)]; }
  double westU(int j, int k) { return m_u[UINDEX(0, j, k, nu, nj, nk)]; }

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

  void divg(Array &g)
  {
    int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
          g[i0] = (m_u[i1 + 1] - m_u[i1]) / (x[i + 1] - x[i]);
          i0++;
          i1++;
        }
        i1++;
      }
    }
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
          g[j0] += (m_v[j1 + 1] - m_v[j1]) / (y[j + 1] - y[j]);
          j0++;
          j1++;
        }
        j1++;
      }
    }
    if (nk > 1) {
      int k0 = 0, k1 = 0;
      for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
          for (int k = 0; k < nk; k++) {
            g[k0] += (m_w[k1 + 1] - m_w[k1]) / (z[k + 1] - z[k]);
            k0++;
            k1++;
          }
          k1++;
        }
      }
    }
  }

  void divg(ArrayU &u, ArrayV &v, Array &g)
  {
    int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    for (int j = 0; j < nj; j++) {
      for (int i = 0; i < ni; i++) {
        g[i0] = (u[i1 + 1] - u[i1]) / (x[i + 1] - x[i]);
        i0++;
        i1++;
      }
      i1++;
    }
    for (int i = 0; i < ni; i++) {
      for (int j = 0; j < nj; j++) {
        g[j0] += (v[j1 + 1] - v[j1]) / (y[j + 1] - y[j]);
        j0++;
        j1++;
      }
      j1++;
    }

  }


  void divg(ArrayU &u, ArrayV &v, ArrayW &w, Array &g)
  {
    int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
          g[i0] = (u[i1 + 1] - u[i1]) / (x[i + 1] - x[i]);
          i0++;
          i1++;
        }
        i1++;
      }
    }
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
          g[j0] += (v[j1 + 1] - v[j1]) / (y[j + 1] - y[j]);
          j0++;
          j1++;
        }
        j1++;
      }
    }
    if (nk > 1) {
      int k0 = 0, k1 = 0;
      for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
          for (int k = 0; k < nk; k++) {
            g[k0] += (w[k1 + 1] - w[k1]) / (z[k + 1] - z[k]);
            k0++;
            k1++;
          }
          k1++;
        }
      }
    }
  }

  void dudx(Array &g)
  {
    int i0 = 0, i1 = 0;
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
          //std::cout << i << " " << j << " " << i0 << " " << i1 << std::endl;
          g[i0] = (u[i1 + 1] - u[i1]) / (x[i + 1] - x[i]);
          i0++;
          i1++;
        }
        i1++;
      }
    }
  }

  void dvdy(Array &g)
  {
    int j0 = 0, j1 = 0;
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
          g[j0] = (v[j1 + 1] - v[j1]) / (y[j + 1] - y[j]);
          j0++;
          j1++;
        }
        j1++;
      }
    }
  }

  void setU(double(*f)(double x, double y))
  {
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        double yy = ym[j];
        for (int i = 0; i < nu; i++) {
          u[UINDEX(i, j, k, nu, nj, nk)] = f(x[i], yy);
        }
      }
    }
  }

  void setV(double(*f)(double x, double y))
  {
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < ni; i++) {
        double xx = xm[i];
        for (int j = 0; j < nv; j++) {
          v[VINDEX(i, j, k, ni, nv, nk)] = f(xx, y[j]);
        }
      }
    }
  }

  void setU(double(*f)(double x, double y, double z))
  {
    if (nk > 1) {
      for (int k = 0; k < nk; k++) {
        double zz = 0.0;
        for (int j = 0; j < nj; j++) {
          double yy = 0.5*(y[j] + y[j + 1]);
          for (int i = 0; i < nu; i++) {
            u[UINDEX(i, j, k, nu, nj, nk)] = f(x[i], yy, zz);
          }
        }
      }
    }
    else {
      for (int j = 0; j < nj; j++) {
        double yy = 0.5*(y[j] + y[j + 1]);
        for (int i = 0; i < nu; i++) {
          u[UINDEX(i, j, 0, nu, nj, nk)] = f(x[i], yy, 0.0);
        }
      }
    }
  }

  void setV(double(*f)(double x, double y, double z))
  {

  }

  void setW(double(*f)(double x, double y, double z))
  {

  }

  void setEastU(double(*f)(double y))
  {
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        //std::cout << 0 << ' ' << j << ' ' << k << ' ' << INDEX(0, j, k, nu, nj, nk) << ' ' << f(ym[j]) << std::endl;
        m_u[INDEX(0, j, k, nu, nj, nk)] = f(ym[j]);
      }
    }
  }

  void setWestU(double(*f)(double y))
  {
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        m_u[INDEX(ni, j, k, nu, nj, nk)] = f(ym[j]);
      }
    }
  }

  void setNorthU(double(*f)(double x))
  {
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < nu; i++) {
        m_u_n[IKINDEX(i, k, nu, nk)] = f(x[i]);
      }
    }
  }

  void setSouthU(double(*f)(double x))
  {
    for (int k = 0; k < nk; k++) {
      for (int i = 0; i < nu; i++) {
        m_u_s[IKINDEX(i, k, nu, nk)] = f(x[i]);
      }
    }
  }

  Array array() {
    return Array(static_cast<StaggeredGrid*>(this));
  }

  ArrayU uArray() {
    return ArrayU(static_cast<StaggeredGrid*>(this));
  }

  ArrayV vArray() {
    return ArrayV(static_cast<StaggeredGrid*>(this));
  }

  ArrayW wArray() {
    return ArrayW(static_cast<StaggeredGrid*>(this));
  }

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

  ArrayU u;
  ArrayV v;
  ArrayW w;
  Array p;

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