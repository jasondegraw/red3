// Copyright (C) 2015-2016 Jason W. DeGraw
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
#include "util.hpp"
#include "staggeredgrid.hpp"
#include "grid.hpp"
#include "arrayops.hpp"
#include <iostream>

namespace red3
{

StaggeredGrid::StaggeredGrid(int ni, int nj, int nk, bool xperi)
  : ni(ni), nj(nj), nk(nk), nu(ni + 1), nv(nj + 1), nw(nk == 1 ? 0 : nk + 1), xperi(xperi)
{
  // Check the inputs
  if(ni < 2) {
    fatal("ni must be greater than 1");
  }
  if(nj < 2) {
    fatal("nj must be greater than 1");
  }

  m_ncells = ni*nj*nk;

  int nuvw = ni+nj+2;
  if(nk > 1) {
    m_nu = nu*nj*nk;
    m_nv = ni*nv*nk;
    m_nw = ni*nj*nw;
    m_w = std::make_unique<double[]>(m_nw);
    //m_w = (double*)callocate(m_nw, sizeof(double), "w velocity");
    //m_z = (double*)callocate(nw, sizeof(double), "z grid");
    //m_zm = (double*)callocate(nz, sizeof(double), "z mid-grid");
  } else {
    m_nu = nu*nj;
    m_nv = ni*nv;
    m_nw = 0;
    //m_w = nullptr;
    //m_z = nullptr;
    //m_zm = nullptr;
  }
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
  if(nw) {
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
  /*
  double dx = 1.0 / (double)nx;
  x.push_back(0.0);
  for(unsigned i = 0; i < nx - 1; i++){
    x.push_back(x[i] + dx);
    xm.push_back(x[i] + 0.5*dx);
  }
  x.push_back(1.0);
  xm.push_back(1.0 - 0.5*dx);
  */

  //s->x = s->y = s->z = 0;
  // s->dx = s->dy = s->dz = 0;
  // s->dx = (real_t*)calloc(s->nx,sizeof(real_t));
  // if(s->dx == NULL)return 0;
  // s->dy = (real_t*)calloc(s->ny,sizeof(real_t));
  // if(s->dy == NULL)return 0;
  // if(nz == 1)  /* Two Dimensions */
    // {
      // s->dim = 2;
      // s->nw  = 0;
     // s->z   = 0;
      // s->w   = 0;
      //s->iv  = (nx+1)*(ny+2);
      //s->iw  = 0;
      //s->nsa = s->iv+(nx+2)*(ny+1);
      // s->nuvw = s->nsa;
      // s->u   = (real_t*)calloc(s->nsa,sizeof(real_t));
      // if(s->u == NULL)return 0;
      // s->v   = s->u + s->iv;
    // }
  // else  /* Three Dimensions */
    // {
      // /* Do it! */
      // /*
          // s->w = (real_t*)malloc(sizeof(real_t)*((nx+2)*(ny+1)*nz));
          // if(s->w == NULL)return 0;
      // */
    // }
  // /*
    // Allocate the pressure array
  // */
  // s->ncell = nx*ny*nz;
  // /* if(s->xperi)s->ncell = (nx+1)*ny*nz; */
  // s->p = (real_t*)calloc(s->ncell,sizeof(real_t));
  // if(s->p == NULL)return 0;
}

double *StaggeredGrid::allocateVariable()
{
  return (double*)callocate(m_ncells, sizeof(double), "cell-centered variable");
}

void StaggeredGrid::divg(Array<StaggeredGrid> &g)
{
  int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      for(int i = 0; i < ni; i++) { 
        g[i0] = (m_u[i1 + 1] - m_u[i1]) / (x[i + 1] - x[i]);
        i0++;
        i1++;
      }
      i1++;
    }
  }
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < ni; i++) {
      for(int j = 0; j < nj; j++) {
        g[j0] += (m_v[j1 + 1] - m_v[j1]) / (y[j + 1] - y[j]);
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
          g[k0] += (m_w[k1 + 1] - m_w[k1]) / (z[k + 1] - z[k]);
          k0++;
          k1++;
        }
        k1++;
      }
    }
  }
}

void StaggeredGrid::divg(ArrayU<StaggeredGrid> &u, ArrayV<StaggeredGrid> &v, Array<StaggeredGrid> &g)
{
  int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
  for(int j = 0; j < nj; j++) {
    for(int i = 0; i < ni; i++) {
      g[i0] = (u[i1 + 1] - u[i1]) / (x[i + 1] - x[i]);
      i0++;
      i1++;
    }
    i1++;
  }
  for(int i = 0; i < ni; i++) {
    for(int j = 0; j < nj; j++) {
      g[j0] += (v[j1 + 1] - v[j1]) / (y[j + 1] - y[j]);
      j0++;
      j1++;
    }
    j1++;
  }

}


void StaggeredGrid::divg(ArrayU<StaggeredGrid> &u, ArrayV<StaggeredGrid> &v, ArrayW<StaggeredGrid> &w, Array<StaggeredGrid> &g)
{
  int i0 = 0, i1 = 0, j0 = 0, j1 = 0;
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      for(int i = 0; i < ni; i++) {
        g[i0] = (u[i1 + 1] - u[i1]) / (x[i + 1] - x[i]);
        i0++;
        i1++;
      }
      i1++;
    }
  }
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < ni; i++) {
      for(int j = 0; j < nj; j++) {
        g[j0] += (v[j1 + 1] - v[j1]) / (y[j + 1] - y[j]);
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
          g[k0] += (w[k1 + 1] - w[k1]) / (z[k + 1] - z[k]);
          k0++;
          k1++;
        }
        k1++;
      }
    }
  }
}

void StaggeredGrid::dudx(Array<StaggeredGrid> &g)
{
  int i0 = 0, i1 = 0;
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      for(int i = 0; i < ni; i++) {
        //std::cout << i << " " << j << " " << i0 << " " << i1 << std::endl;
        g[i0] = (m_u[i1 + 1] - m_u[i1]) / (x[i + 1] - x[i]);
        i0++;
        i1++;
      }
      i1++;
    }
  }
}

void StaggeredGrid::dvdy(Array<StaggeredGrid> &g)
{
  int j0 = 0, j1 = 0;
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < ni; i++) {
      for(int j = 0; j < nj; j++) {
        g[j0] = (m_v[j1 + 1] - m_v[j1]) / (y[j + 1] - y[j]);
        j0++;
        j1++;
      }
      j1++;
    }
  }
}

//void StaggeredGrid::setU(std::function<double(double, double)> f)
void StaggeredGrid::setU(double(*f)(double x, double y))
{
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      double yy = ym[j];
      for(int i = 0; i < nu; i++) {
        m_u[UINDEX(i, j, 0, nu, nj, nk)] = f(x[i], yy);
      }
    }
  }
}

void StaggeredGrid::setV(double(*f)(double x, double y))
{
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < ni; i++) {
      double xx = xm[i];
      for(int j = 0; j < nv; j++) {
        m_v[VINDEX(i, j, k, ni, nv, nk)] = f(xx, y[j]);
      }
    }
  }
}

//void StaggeredGrid::setU(std::function<double(double, double, double)> f)
void StaggeredGrid::setU(double(*f)(double x, double y, double z))
{
  if(nw) {
    for(int k = 0; k < nk; k++) {
      double zz = 0.0;
      for(int j = 0; j < nj; j++) {
        double yy = 0.5*(y[j] + y[j + 1]);
        for(int i = 0; i < nu; i++) {
          m_u[UINDEX(i, j, k, nu, nj, nk)] = f(x[i], yy, zz);
        }
      }
    }
  } else {
    for(int j = 0; j < nj; j++) {
      double yy = 0.5*(y[j] + y[j + 1]);
      for(int i = 0; i < nu; i++) {
        m_u[UINDEX(i, j, 0, nu, nj, nk)] = f(x[i], yy, 0.0);
      }
    }
  }
}

void StaggeredGrid::setV(double(*f)(double x, double y, double z))
{

}

void StaggeredGrid::setW(double(*f)(double x, double y, double z))
{

}

void StaggeredGrid::setEastU(double(*f)(double y))
{
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      m_u[UINDEX(0, j, k, ni, nv, nk)] = f(ym[j]);
    }
  }
}

void StaggeredGrid::setWestU(double(*f)(double y))
{
  for(int k = 0; k < nk; k++) {
    for(int j = 0; j < nj; j++) {
      m_u[UINDEX(nu, j, k, ni, nv, nk)] = f(ym[j]);
    }
  }
}

void StaggeredGrid::setNorthU(double(*f)(double x))
{
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < nu; i++) {
      m_u_n[IKINDEX(i, k, nu, nk)] = f(x[i]);
    }
  }
}

void StaggeredGrid::setSouthU(double(*f)(double x))
{
  for(int k = 0; k < nk; k++) {
    for(int i = 0; i < nu; i++) {
      m_u_s[IKINDEX(i, k, nu, nk)] = f(x[i]);
    }
  }
}

}
