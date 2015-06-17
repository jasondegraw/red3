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
#include "util.h"
#include "staggeredmodel.h"

namespace red3
{

StaggeredModel::StaggeredModel(unsigned nx, unsigned ny, unsigned nz)
{
  int i;
  // Check the inputs
  if(nx <= 1) {
    fatal("nx must be greater than 1");
  }
  if(ny <= 1) {
    fatal("ny must be greater than 1");
  }
  if(nz < 1) {
    fatal("nz must be greater than 0");
  }

  m_nx = nx;
  m_ny = ny;
  m_nz = nz;
  m_nu = nx+1;
  m_nv = ny+1;
  m_nw = nz+1;
  unsigned nuvw = nx+ny+2;
  if(nz > 1) {
    nuvw += nz+1;
  }

}

}
