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
#ifndef RED3_STAGGEREDGRID_HPP
#define RED3_STAGGEREDGRID_HPP
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

//using index_t = int;

class RED3_API StaggeredGrid
{
public:
  using ArrayM = ChildArrayUVW<StaggeredGrid>;
  using ArrayP = ChildArray<StaggeredGrid>;
  using ArrayU = ChildArrayU<StaggeredGrid>;
  using ArrayV = ChildArrayV<StaggeredGrid>;
  using ArrayW = ChildArrayW<StaggeredGrid>;

  //using index_t = int;

  StaggeredGrid(Grid1D& x, Grid1D& y, Grid1D& z = Grid1D::one(), bool xperi = false);

  virtual ~StaggeredGrid()
  {}

  double northU(int i, int k) { return m_u_n[IKINDEX(i, k, nu, nk)]; }
  double southU(int i, int k) { return m_u_s[IKINDEX(i, k, nu, nk)]; }

  double eastU(int j, int k) { return u[UINDEX(ni, j, k, nu, nj, nk)]; }
  double westU(int j, int k) { return u[UINDEX(0, j, k, nu, nj, nk)]; }

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

  void divg(ArrayP &g)
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

  void divg(ArrayU &u, ArrayV &v, ArrayP &g)
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


  void divg(ArrayU &u, ArrayV &v, ArrayW &w, ArrayP &g)
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

  void dudx(ArrayP &g)
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

  void dvdy(ArrayP &g)
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
        double yy = y.m(j);
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
        double xx = x.m(i);
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
        u[INDEX(ni, j, k, nu, nj, nk)] = f(y.m(j));
      }
    }
  }

  void setWestU(double(*f)(double y))
  {
    for (int k = 0; k < nk; k++) {
      for (int j = 0; j < nj; j++) {
        u[INDEX(0, j, k, nu, nj, nk)] = f(y.m(j));
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

  ArrayM maxArray() {
    return ArrayM(static_cast<StaggeredGrid*>(this));
  }

  ArrayP pArray() {
    return ArrayP(static_cast<StaggeredGrid*>(this));
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

  void divg(double *g);

  const index_t ni;
  const index_t nj;
  const index_t nk;

  const index_t nu;
  const index_t nv;
  const index_t nw;

  const index_t nnu;
  const index_t nnv;
  const index_t nnw;
  const index_t ncells;

  const bool xperi;
  const bool two_dimensional;

  const Grid1D x;
  const Grid1D y;
  const Grid1D z;

  //std::vector<double> x;
  //std::vector<double> dx;
  //std::vector<double> y;
  //std::vector<double> dy;
  //std::vector<double> z;
  //std::vector<double> dz;
  //std::vector<double> xm;
  //std::vector<double> ym;
  //std::vector<double> zm;

  ArrayU u;
  ArrayV v;
  ArrayW w;
  ArrayP p;

  double bfx, bfy, bfz;

protected:

  int m_dimension;

  double *m_u_n;
  double *m_u_s;

};

}

#endif