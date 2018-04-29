// Copyright (C) 2015-2017 Jason W. DeGraw
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
#include <stdlib.h>
#include "upwind.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
  red3::upwind::StaggeredIncompressibleSteadyFlow solver(100.0, 1.0, 9, 9);
  red3::Array<red3::upwind::StaggeredIncompressibleSteadyFlow> pstar = solver.array();
  red3::ArrayU<red3::upwind::StaggeredIncompressibleSteadyFlow> ustar = solver.uArray();
  red3::ArrayV<red3::upwind::StaggeredIncompressibleSteadyFlow> vstar = solver.vArray();

  solver.setEastU([](double x) {return 4*x*(1.0 - x); });
  solver.setWestU([](double x) {return 4*x*(1.0 - x); });
  for (int j = 0; j < solver.nj; j++) {
    std::cout << solver.xm[j] << ' ' << solver.u(0, j, 0) << ' ' << solver.u(solver.ni, j, 0) << std::endl;
  }

  std::cout << solver.ni << ' ' << solver.nj << ' ' << solver.nk << std::endl;
  std::cout << solver.nu << ' ' << solver.nv << ' ' << solver.nw << std::endl;
  std::cout << solver.ae.size() << ' ' << solver.aw.size() << ' ' << solver.an.size() << ' ' << solver.as.size() << std::endl;
  std::cout << solver.ae.parent() << ' ' << &solver << std::endl;

  solver.setupU();

  for (int j = 0; j < solver.nj; j++) {
    std::cout << j << ' ' << solver.xm[j] << ' ' << solver.eastU(j, 0) << std::endl; // ' ' << solver.u(solver.ni, j, 0) << std::endl;
  }

  //auto us = solver.uArray();

  return EXIT_SUCCESS;
}
