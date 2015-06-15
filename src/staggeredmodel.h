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
private:
  double m_dt;
  double *x;
  double *xm;
  double *y;
  double *ym;
  double *z;
  double *zm;
  double *u;
  double *v;
  double *w;
};

}
