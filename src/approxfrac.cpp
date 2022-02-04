// Copyright (C) 2015-2022 Jason W. DeGraw
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
#include "red3/approxfrac.hpp"
#include "red3/grid.hpp"
#include "red3/arrayops.hpp"
#include "red3/utilities.hpp"
#include <iostream>

namespace red3 {
namespace approxfrac {

ViscousOperatorU::ViscousOperatorU(StaggeredGrid *parent, double dt, double Re) : StaggeredGrid::ArrayU(parent) //, last(parent)
{
  double mul = dt / Re;
  m_ax.resize(parent->nu-2);
  m_bx.resize(parent->nu-2);
  m_cx.resize(parent->nu-2);
  if (parent->uniform_x) {
    double dx = parent->dx[0];
    double rdx2 = mul / (dx*dx);
    for (index_t i = 0; i < parent->nu-2; ++i) {
      m_ax[i] = m_cx[i] = -rdx2;
      m_bx[i] = 2.0*rdx2;
    }
  } else {
    // TODO
  }

  m_ay.resize(parent->nj);
  m_by.resize(parent->nj);
  m_cy.resize(parent->nj);
  if (parent->uniform_y) {
    double dy = parent->dy[0];
    double rdy2 = mul / (dy*dy);
    for (index_t i = 0; i < parent->nj; ++i) {
      m_ay[i] = m_cy[i] = -rdy2;
      m_bx[i] = 2.0*rdy2;
    }
  } else {
    // TODO
  }

  if (!parent->two_dimensional) {
    // TODO
  }

}

void ViscousOperatorU::compute(const StaggeredGrid::ArrayU &u)
{

}


ViscousOperatorV::ViscousOperatorV(StaggeredGrid *parent, double dt, double Re) : StaggeredGrid::ArrayV(parent) //, last(parent)
{
  double mul = dt / Re;
  m_ax.resize(parent->ni);
  m_bx.resize(parent->ni);
  m_cx.resize(parent->ni);
  if (parent->uniform_x) {
    double dx = parent->dx[0];
    double rdx2 = mul / (dx*dx);
    for (index_t i = 0; i < parent->ni; ++i) {
      m_ax[i] = m_cx[i] = -rdx2;
      m_bx[i] = 2.0*rdx2;
    }
  } else {
    // TODO
  }

  m_ay.resize(parent->nj);
  m_by.resize(parent->nj);
  m_cy.resize(parent->nj);
  if (parent->uniform_y) {
    double dy = parent->dy[0];
    double rdy2 = mul / (dy*dy);
    for (index_t i = 0; i < parent->nj; ++i) {
      m_ay[i] = m_cy[i] = -rdy2;
      m_bx[i] = 2.0*rdy2;
    }
  } else {
    // TODO
  }

  if (!parent->two_dimensional) {
    // TODO
  }

}

void ViscousOperatorV::compute(const StaggeredGrid::ArrayV &u)
{

}

void IsothermalFlow::convec(ArrayM &h1, ArrayM &h2, ArrayM &h3, index_t i1, index_t i2, index_t j1, index_t j2, 
  index_t k1, index_t k2)
{
  //
  //  *******************************************************************
  //
  //     Compute the convective term
  //
  //     The original version was correct only for uniform grids.This
  //     version is correct.
  //
  //     Author : Jason W.DeGraw(06-22-2005)
  //  *******************************************************************
  //   integer nx, ny, nz, i1, i2, j1, j2, k1, k2
  //  real * 8 u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz),
  //   $     h1(nx, ny, nz), h2(nx, ny, nz), h3(nx, ny, nz)
  //   real * 8 y(ny), ym(ny)
  //
  //   integer imax, jmax, kmax, imaxm1, imaxp1, jmaxm1, jmaxp1,
  //   $     nxfft, nzfft
  // common / data9 / imax, jmax, kmax, imaxm1, imaxp1, jmaxm1, jmaxp1,
  //   $     nxfft, nzfft
  //   real * 8 delta1, delta3, re
  // common / mesh / delta1, delta3, re
  //   real * 8 bfx, bfy, bfz
  //   common / bforce / bfx, bfy, bfz
  //   c
  //   real * 8 ue2x4, uw2x4, un, us, vnx2, vsx2,
  //   $     ue, uw, vex2, vwx2, vn2x4, vs2x4, vbx2, vfx2, wex2, wwx2,
  //   $     ufx2, ubx2, wfx2, wbx2, uex2, uwx2, wn, ws, wf2x4, wb2x4
  //   c
  //   real * 8 cz, cx, cxv, cy, cym, cyn, cys
  //   integer i, j, k, km1, kp1, im1, ip1, jm1, jp1
  //
  double cx = 1.0; // 0.25 / delta1;
  double cxv = 1.0; // 0.5 / delta1;
  double cz = 1.0; // 0.25 / delta3;

  if (two_dimensional) {
    //
    //.... 2D version
    //
    index_t k = 1;
    for (index_t j = j1; j <= j2; ++j) {
      index_t jm1 = j - 1;
      index_t jp1 = j + 1;
      double cy = 0.5 / (y[j] - y[jm1]);
      double cym = 0.25 / (ym[jp1] - ym[j]);
      double cyn = (y[j] - ym[j]) / (ym[jp1] - ym[j]);
      double cys = (y[jm1] - ym[jm1]) / (ym[j] - ym[jm1]);
      for (index_t i = i1; i <= i2; ++i) {
        //do i = i1, i2          !2, imax
        index_t im1 = i - 1;
        index_t ip1 = i + 1;
        double ue2x4 = power2(u(ip1, j, k) + u(i, j, k));
        double uw2x4 = power2(u(im1, j, k) + u(i, j, k));
        double un = cyn*u(i, jp1, k) + (1.0 - cyn)*u(i, j, k);
        double us = cys*u(i, j, k) + (1.0 - cys)*u(i, jm1, k);
        double vnx2 = (v(ip1, j, k) + v(i, j, k));
        double vsx2 = (v(ip1, jm1, k) + v(i, jm1, k));
        //
        //....x - momentum equation.
        //
        h1(i, j, k) = bfx - cx*(ue2x4 - uw2x4) - cy*(un*vnx2 - us*vsx2);
        //
        double ue = un;
        double uw = cyn*u(im1, jp1, k) + (1.0 - cyn)*u(im1, j, k);
        double vex2 = vnx2;
        double vwx2 = (v(i, j, k) + v(im1, j, k));
        double vn2x4 = power2(v(i, jp1, k) + v(i, j, k));
        double vs2x4 = power2(v(i, jm1, k) + v(i, j, k));
        //
        //....y - momentum equation.
        //
        h2(i, j, k) = bfy - cxv*(ue*vex2 - uw*vwx2) - cym*(vn2x4 - vs2x4);
      }
    }
  } else {
    //
    //.... 3D version
    //
    for (index_t k = k1; k <= k2; ++k) {
      //do k = k1, k2             !1, kmax
      index_t km1 = k - 1;
      index_t kp1 = k + 1;
      if (k == k1) {
        km1 = k2;
      } else if (k == k2) {
        kp1 = k1;
      }
      for (index_t j = j1; j <= j2; ++j) {
        //do j = j1, j2          !2, jmax
        index_t jm1 = j - 1;
        index_t jp1 = j + 1;
        double cy = 0.5 / (y[j] - y[j - 1]);
        double cym = 0.25 / (ym[jp1] - ym[j]);
        double cyn = (y[j] - ym[j]) / (ym[jp1] - ym[j]);
        double cys = (y[jm1] - ym[jm1]) / (ym[j] - ym[jm1]);
        for (index_t i = i1; i <= i2; ++i) {
          //do i = i1, i2       !2, imax
          index_t im1 = i - 1;
          index_t ip1 = i + 1;
          double ue2x4 = power2(u(ip1, j, k) + u(i, j, k));
          double uw2x4 = power2(u(im1, j, k) + u(i, j, k));
          double un = cyn*u(i, jp1, k) + (1.0 - cyn)*u(i, j, k);
          double us = cys*u(i, j, k) + (1.0 - cys)*u(i, jm1, k);
          double vnx2 = (v(ip1, j, k) + v(i, j, k));
          double vsx2 = (v(ip1, jm1, k) + v(i, jm1, k));
          double ufx2 = (u(i, j, k) + u(i, j, kp1));
          double ubx2 = (u(i, j, k) + u(i, j, km1));
          double wfx2 = (w(i, j, k) + w(ip1, j, k));
          double wbx2 = (w(i, j, k) + w(ip1, j, km1));
          //
          //....x - momentum equation.
          //
          h1(i, j, k) = bfx - cx*(ue2x4 - uw2x4) - cy*(un*vnx2 - us*vsx2) - cz*(ufx2*wfx2 - ubx2*wbx2);
          //
          double ue = un;
          double uw = cyn*u(im1, jp1, k) + (1.0 - cyn)*u(im1, j, k);
          double vex2 = vnx2;
          double vwx2 = (v(i, j, k) + v(im1, j, k));
          double vn2x4 = power2(v(i, jp1, k) + v(i, j, k));
          double vs2x4 = power2(v(i, jm1, k) + v(i, j, k));
          double vfx2 = (v(i, j, k) + v(i, j, kp1));
          double vbx2 = (v(i, j, k) + v(i, j, km1));
          wfx2 = (w(i, j, k) + w(i, jp1, k));
          wbx2 = (w(i, j, km1) + w(i, jp1, km1));
          //
          //....y - momentum equation.
          //
          h2(i, j, k) = bfy - cxv*(ue*vex2 - uw*vwx2) - cym*(vn2x4 - vs2x4) - cz*(vfx2*wfx2 - vbx2*wbx2);
          //
          double uex2 = ufx2;
          double uwx2 = (u(im1, j, k) + u(im1, j, kp1));
          double wex2 = wfx2;
          double wwx2 = (w(i, j, k) + w(im1, j, k));
          vnx2 = vfx2;
          vsx2 = (v(i, jm1, k) + v(i, jm1, kp1));
          double wn = cyn*w(i, jp1, k) + (1.0 - cyn)*w(i, j, k);
          double ws = cys*w(i, j, k) + (1.0 - cys)*w(i, jm1, k);
          double wf2x4 = power2(w(i, j, kp1) + w(i, j, k));
          double wb2x4 = power2(w(i, j, km1) + w(i, j, k));
          //
          //....z - momentum equation
          //
          h3(i, j, k) = bfz - cx*(uex2*wex2 - uwx2*wwx2) - cy*(vnx2*wn - vsx2*ws) - cz*(wf2x4 - wb2x4);
        }
      }
    }
  }
}

/*
void IsothermalFlow::visc(u, g, ax, bx, cx, ay, by, cy, i1, i2, j1, j2,
  $     ncase, nx, ny, nz)
{
  //  *******************************************************************
  //     This is basically the original version, but it has been reworked
  //     a bit to get it to work better.In particular :
  //       - the pressure gradient addition has been moved elsewhere
  //       - the arguments reordered and some added
  //       - the damping stuff gone
  //
  //     OK, so it is not very much like the original.
  //
  //     "Author": Jason W.DeGraw(09 - 16 - 2005)
  //
  //  *******************************************************************
  // implicit none !double precision(a - h, o - z)
  // integer i1, i2, j1, j2, ncase, nx, ny, nz
  // real * 8 u(nx, ny, nz), g(nx, ny, nz)
  // real * 8 ax(nx), bx(nx), cx(nx)
  // real * 8 ay(ny), by(ny), cy(ny)
  //
  // include 'common.f'
  // integer imax, jmax, kmax, imaxm1, imaxp1, jmaxm1, jmaxp1,
  // &             nxfft, nzfft, nend, ntime, jmid
  // common / data9 / imax, jmax, kmax, imaxm1, imaxp1, jmaxm1, jmaxp1,
  // &             nxfft, nzfft
  // double precision delta1, delta3, re, dt
  // common / step / dt, nend, ntime
  // common / mesh / delta1, delta3, re, jmid
  // integer icpr, iback, nread, nsave, nbc, iturb
  // common / flag / icpr, iback, nread, nsave, nbc, iturb
  //
  // real * 8 amx(nnx), bmx(nnx), cmx(nnx),
  // &amz(nnx), bmz(nnx), cmz(nnx), amy(nnx), bmy(nnx), cmy(nnx),
  // &f(nnx), q(nnx), s(nnx)
  // common / room / amx, bmx, cmx, amz, bmz, cmz, amy, bmy, cmy, f, q, s
  //
  // integer i, j, k, km1, kp1, im1, ip1
  // real * 8 xmul, ymul, az

  double xmul = dt / re;
  double ymul = dt / re;
  double az = dt / (re*delta3*delta3);
  double xmul = dt / re;
c
do i = i1, i2
amx(i) = ax(i)*xmul
bmx(i) = bx(i)*xmul
cmx(i) = cx(i)*xmul
enddo
do i = 1, kmax
amz(i) = az
bmz(i) = -2.0d0*az
cmz(i) = az
enddo
do i = j1, j2
amy(i) = ay(i)*ymul
bmy(i) = by(i)*ymul
cmy(i) = cy(i)*ymul
enddo
c
if (nbc.eq. 7)then
do k = 1, kmax
km1 = k - 1
kp1 = k + 1
if (k.eq.1) km1 = kmax
if (k.eq.kmax) kp1 = 1
i = i1
ip1 = i + 1
im1 = i2
do j = j1, j2
g(i, j, k) =
$               amx(i)*u(im1, j, k) + bmx(i)*u(i, j, k) + cmx(i)*u(ip1, j, k)
$ + amy(j)*u(i, j - 1, k) + bmy(j)*u(i, j, k) + cmy(j)*u(i, j + 1, k)
$ + amz(k)*u(i, j, km1) + bmz(k)*u(i, j, k) + cmz(k)*u(i, j, kp1)

enddo
do i = i1 + 1, i2 - 1
im1 = i - 1
ip1 = i + 1
do j = j1, j2
g(i, j, k) =
$               amx(i)*u(im1, j, k) + bmx(i)*u(i, j, k) + cmx(i)*u(ip1, j, k)
$ + amy(j)*u(i, j - 1, k) + bmy(j)*u(i, j, k) + cmy(j)*u(i, j + 1, k)
$ + amz(k)*u(i, j, km1) + bmz(k)*u(i, j, k) + cmz(k)*u(i, j, kp1)

enddo
enddo
i = i2
ip1 = i1
im1 = i - 1
do j = j1, j2
g(i, j, k) =
$               amx(i)*u(im1, j, k) + bmx(i)*u(i, j, k) + cmx(i)*u(ip1, j, k)
$ + amy(j)*u(i, j - 1, k) + bmy(j)*u(i, j, k) + cmy(j)*u(i, j + 1, k)
$ + amz(k)*u(i, j, km1) + bmz(k)*u(i, j, k) + cmz(k)*u(i, j, kp1)

enddo
enddo
else
do k = 1, kmax
km1 = k - 1
kp1 = k + 1
if (k.eq.1) km1 = kmax
if (k.eq.kmax) kp1 = 1
do i = i1, i2
im1 = i - 1
ip1 = i + 1
do j = j1, j2
g(i, j, k) =
$               amx(i)*u(im1, j, k) + bmx(i)*u(i, j, k) + cmx(i)*u(ip1, j, k)
$ + amy(j)*u(i, j - 1, k) + bmy(j)*u(i, j, k) + cmy(j)*u(i, j + 1, k)
$ + amz(k)*u(i, j, km1) + bmz(k)*u(i, j, k) + cmz(k)*u(i, j, kp1)

enddo
enddo
enddo
endif
c
return
end
*/

  /*

StaggeredGrid::StaggeredGrid(int ni, int nj, int nk, bool xperi)
  : ni(ni), nj(nj), nk(nk), nu(ni + 1), nv(nj + 1), nw(nk == 1 ? 1 : nk + 1), xperi(xperi)
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
  if(nk > 1) {
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

  //double dx = 1.0 / (double)nx;
  //x.push_back(0.0);
  //for(unsigned i = 0; i < nx - 1; i++){
  //  x.push_back(x[i] + dx);
  //  xm.push_back(x[i] + 0.5*dx);
  //}
  //x.push_back(1.0);
  //xm.push_back(1.0 - 0.5*dx);

  //s->x = s->y = s->z = 0;
  // s->dx = s->dy = s->dz = 0;
  // s->dx = (real_t*)calloc(s->nx,sizeof(real_t));
  // if(s->dx == NULL)return 0;
  // s->dy = (real_t*)calloc(s->ny,sizeof(real_t));
  // if(s->dy == NULL)return 0;
  // if(nz == 1)  //Two Dimensions
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
  // else  // Three Dimensions
    // {
      // Do it!
      //
          // s->w = (real_t*)malloc(sizeof(real_t)*((nx+2)*(ny+1)*nz));
          // if(s->w == NULL)return 0;
      //
    // }
  //
    // Allocate the pressure array
  //
  // s->ncell = nx*ny*nz;
  // if(s->xperi)s->ncell = (nx+1)*ny*nz;
  // s->p = (real_t*)calloc(s->ncell,sizeof(real_t));
  // if(s->p == NULL)return 0;
}

double *StaggeredGrid::allocateVariable()
{
  return (double*)callocate(m_ncells, sizeof(double), "cell-centered variable");
}





//void StaggeredGrid::setU(std::function<double(double, double)> f)


//void StaggeredGrid::setU(std::function<double(double, double, double)> f)
void StaggeredGrid::setU(double(*f)(double x, double y, double z))
{
  if(nk > 1) {
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

*/
}
}
