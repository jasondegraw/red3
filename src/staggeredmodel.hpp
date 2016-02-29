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

#include "red3api.h"

namespace red3
{

class RED3_API StaggeredModel
{
public:
  StaggeredModel(unsigned nx, unsigned ny, unsigned nz=1);
private:
  unsigned m_nx;
  unsigned m_ny;
  unsigned m_nz;
  unsigned m_nu;
  unsigned m_nv;
  unsigned m_nw;

  double m_reynum;
  double m_dt;
  double *m_x;
  double *m_xm;
  double *m_y;
  double *m_ym;
  double *m_z;
  double *m_zm;
  double *m_u;
  double *m_v;
  double *m_w;
  double *m_p;
};

}
