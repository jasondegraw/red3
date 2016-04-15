// Copyright (C) 2015-2015 Jason W. DeGraw
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

namespace red3
{

StaggeredGrid::StaggeredGrid(unsigned ni, unsigned nj, unsigned nk) : nx(ni), ny(nj), nz(nk), nu(ni + 1), nv(nj + 1), nw(nk==1 ? 0 : nk+1)
{
  int i;
  // Check the inputs
  if(nx <= 1) {
    fatal("nx must be greater than 1");
  }
  if(ny <= 1) {
    fatal("ny must be greater than 1");
  }

  m_ncells = nx*ny*nz;

  unsigned nuvw = nx+ny+2;
  if(nz > 1) {
    m_nu = nu*ny*nz;
    m_nv = nx*nv*nz;
    m_nw = nx*ny*nw;
    m_w = (double*)callocate(m_nw, sizeof(double), "w velocity");
    //m_z = (double*)callocate(nw, sizeof(double), "z grid");
    //m_zm = (double*)callocate(nz, sizeof(double), "z mid-grid");
  } else {
    m_nu = nu*ny;
    m_nv = nx*nv;
    m_nw = 0;
    m_w = nullptr;
    //m_z = nullptr;
    //m_zm = nullptr;
  }
  m_u = (double*)callocate(m_nu, sizeof(double), "u velocity");
  m_v = (double*)callocate(m_nv, sizeof(double), "v velocity");
  //m_x = (double*)callocate(nu, sizeof(double), "x grid");
  //m_y = (double*)callocate(nv, sizeof(double), "y grid");
  //m_xm = (double*)callocate(nx, sizeof(double), "x mid-grid");
  //m_ym = (double*)callocate(ny, sizeof(double), "y mid-grid");
  m_p = (double*)callocate(m_ncells, sizeof(double), "pressure");

  // Make a simple grid
  Uniform xg(nu);
  x = xg.grid();
  xm = xg.midgrid();
  Uniform yg(nv);
  y = yg.grid();
  ym = yg.midgrid();
  if(nw){
    Uniform zg(nv);
    z = zg.grid();
    zm = zg.midgrid();
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

StaggeredGrid::~StaggeredGrid()
{
  if(m_nw) {
    free(m_w);
  }
  free(m_u);
  free(m_v);
  free(m_p);
}

//void StaggeredGrid::setU(std::function<double(double, double)> f)
void StaggeredGrid::setU(double(*f)(double x, double y))
{
  for(unsigned k = 0; k < nz; k++) {
    for(unsigned j = 0; j < ny; j++) {
      double yy = ym[j];
      for(unsigned i = 0; i < nu; i++) {
        m_u[UINDEX(i, j, 0, nu, ny, nz)] = f(x[i], yy);
      }
    }
  }
}

void StaggeredGrid::setV(double(*f)(double x, double y))
{
  for(unsigned k = 0; k < nz; k++) {
    for(unsigned i = 0; i < nx; i++) {
      double xx = xm[i];
      for(unsigned j = 0; j < nv; j++) {
        m_u[UINDEX(i, j, k, nu, ny, nz)] = f(xx, y[i]);
      }
    }
  }
}

//void StaggeredGrid::setU(std::function<double(double, double, double)> f)
void StaggeredGrid::setU(double(*f)(double x, double y, double z))
{
  if(nw) {
    for(unsigned k = 0; k < nz; k++) {
      double zz = 0.0;
      for(unsigned j = 0; j < ny; j++) {
        double yy = 0.5*(y[j] + y[j + 1]);
        for(unsigned i = 0; i < nu; i++) {
          m_u[UINDEX(i, j, k, nu, ny, nz)] = f(x[i], yy, zz);
        }
      }
    }
  } else {
    for(unsigned j = 0; j < ny; j++) {
      double yy = 0.5*(y[j] + y[j+1]);
      for(unsigned i = 0; i < nu; i++) {
        m_u[UINDEX(i, j, 0, nu, ny, nz)] = f(x[i], yy, 0.0);
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

double *StaggeredGrid::allocateVariable()
{
  return (double*)callocate(m_ncells, sizeof(double), "cell-centered variable");
}

void StaggeredGrid::divg(Array<StaggeredGrid> &g)
{
  //Array<StaggeredGrid> g(this);
  unsigned i0 = 0, i1 = 0;
  for(unsigned k = 0; k < nz; k++) {
    for(unsigned j = 0; j < ny; j++) {
      for(unsigned i = 0; i < nx; i++) { 
        g[i0] = (m_u[i1 + 1] - m_u[i1]) / dx[i];
        i0++;
        i1++;
      }
      i1++;
    }
  }
  //return g;
}

void StaggeredGrid::dudx(Array<StaggeredGrid> &g)
{
  //Array<StaggeredGrid> g(this);
  unsigned i0 = 0, i1 = 0;
  for(unsigned k = 0; k < nz; k++) {
    for(unsigned j = 0; j < ny; j++) {
      for(unsigned i = 0; i < nx; i++) {
        g[i0] = (m_u[i1 + 1] - m_u[i1]) / dx[i];
        i0++;
        i1++;
      }
      i1++;
    }
  }
  //return g;
}

}
