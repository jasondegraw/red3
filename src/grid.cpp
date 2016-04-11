// Copyright (C) 2015 Jason W. DeGraw
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

Uniform::Uniform(unsigned n, double L, bool lengthIsDelta)
{
  m_n = n;
  if(lengthIsDelta){
    m_dx = L;
    m_L = L*n;
  } else {
    m_L = L;
    m_dx = L/(double)(n-1);
  }
}

std::vector<double> Uniform::grid()
{
  std::vector<double> g = {0.0};
  for(unsigned i=1; i<m_n-1; i++) {
    g.push_back(i*m_dx);
  }
  g.push_back(m_L);
  return g;
}

std::vector<double> Uniform::midgrid()
{
  std::vector<double> g;
  for(unsigned i=1; i<m_n; i++) {
    g.push_back((i-0.5)*m_dx);
  }
  return g;
}
