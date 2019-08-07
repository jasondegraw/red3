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

#include "staggeredgrid.hpp"
#include "array.hpp"
#include "catch.hpp"

TEST_CASE("Basic 2D Arrays", "[Array]")
{
  red3::StaggeredGrid grid(4, 5);
  red3::StaggeredGrid::ArrayU array(&grid);
  for(auto j = 0; j < grid.nj; j++) {
    for(auto i = 0; i < grid.nu; i++) {
      array(i, j, 0) = i;
    }
  }
  std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0,
    0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0};
  for(auto i = 0; i < 25; i++) {
    REQUIRE(v[i] == array[i]);
    array[i] = (double)i;
  }
  double value = 0.0;
  for(auto i = 0; i < 25; i++) {
    REQUIRE(value == array[i]);
    value += 1.0;
  }
}

