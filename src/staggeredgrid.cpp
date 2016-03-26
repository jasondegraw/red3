// Copyright (C) 2015  Jason W. DeGraw
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

namespace red3
{

StaggeredGrid::StaggeredGrid(unsigned nx, unsigned ny, unsigned nz)
{
  int i;
  // Check the inputs
  if(nx <= 1) {
    fatal("nx must be greater than 1");
  }
  if(ny <= 1) {
    fatal("ny must be greater than 1");
  }

  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  nu = nx+1;
  nv = ny+1;
  nw = nz+1;

  unsigned nuvw = nx+ny+2;
  if(nz > 1) {
    m_nu = nu*ny*nz;
    m_nv = nx*nv*nz;
    m_nw = nx*ny*nw;
    m_w = (double*)calloc(m_nw, sizeof(double));
  } else {
    m_nu = nu*ny;
    m_nv = nx*nv;
    m_nw = 0;
    m_w = nullptr;
  }
  m_u = (double*)calloc(m_nu, sizeof(double));
  m_v = (double*)calloc(m_nv, sizeof(double));
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
}

}
