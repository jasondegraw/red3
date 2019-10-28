// Copyright (C) 2015-2018 Jason W. DeGraw
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

#include "upwind.hpp"
#include "catch.hpp"
#include <iostream>
#include <string>

TEST_CASE("Basic 2D Upwinding", "[upwind]")
{
  red3::upwind::StaggeredIncompressibleSteadyFlow grid(100.0, 1.0, red3::Grid1D(4), red3::Grid1D(5));
  REQUIRE(4 == grid.ni);
  REQUIRE(5 == grid.nu);
  REQUIRE(5 == grid.nj);
  REQUIRE(6 == grid.nv);
  REQUIRE(1 == grid.nk);
  REQUIRE(1 == grid.nw);
  // Check grid
  REQUIRE(5 == grid.x.size());
  std::vector<double> x = { 0.0, 0.25, 0.5, 0.75, 1.0 };
  for(auto i = 0; i < grid.nu; i++) {
    REQUIRE(x[i] == grid.x[i]);
  }
  REQUIRE(4 == grid.x.midgrid().size());
  std::vector<double> xm = { 0.125, 0.375, 0.625, 0.875 };
  for(auto i = 0; i < grid.ni; i++) {
    REQUIRE(xm[i] == grid.x.m(i));
  }
  REQUIRE(6 == grid.y.size());
  std::vector<double> y = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
  for(auto i = 0; i < grid.nv; i++) {
    REQUIRE(y[i] == Approx(grid.y[i]));
  }
  REQUIRE(5 == grid.y.midgrid().size());
  std::vector<double> ym = { 0.1, 0.3, 0.5, 0.7, 0.9 };
  for(auto i = 0; i < grid.ni; i++) {
    REQUIRE(ym[i] == Approx(grid.y.m(i)));
  }
  // BC Setters
  grid.setNorthU([](double x){return 1.0;});
  for(auto i = 0; i < grid.nu; i++) {
    REQUIRE(1.0 == grid.northU(i, 0)); // << ("Index: " + std::to_string(i));
  }

}

TEST_CASE("1D Tests", "[upwind]")
{
  red3::upwind::StaggeredIncompressibleSteadyFlow grid(100.0, 1.0, red3::Grid1D(5), red3::Grid1D(2));
  REQUIRE(5 == grid.ni);
  REQUIRE(6 == grid.nu);
  REQUIRE(2 == grid.nj);
  REQUIRE(3 == grid.nv);
  REQUIRE(1 == grid.nk);
  REQUIRE(1 == grid.nw);
}

