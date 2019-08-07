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
#include"grid.hpp"

namespace red3 {

Grid1D::Grid1D(index_t n) : generator(Generator::UniformN), N(n), uniform(true)
{
  double d{ 1.0 / static_cast<double>(N - 1) };
  m_dx.resize(N - 1);
  m_x.resize(N);
  m_xm.resize(N - 1);
  m_x[0] = 0.0;
  m_xm[0] = 0.5*d;
  m_dx[0] = 0.0;
  for (index_t i = 1; i < N - 1; ++i) {
    m_dx[i] = d;
    m_x[i] = m_x[i - 1] + d;
    m_xm[i] = m_xm[i - 1] + d;
  }
  m_x[N - 1] = 1.0;
}

}