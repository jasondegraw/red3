// Copyright (C) 2015-2019 Jason W. DeGraw
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
#include "catch.hpp"

#include "grid.hpp"

TEST_CASE("Uniform Grid 1D", "[grid]")
{
  red3::Uniform1D uniform_x(0.25, 4);
  std::optional<red3::Grid1D> opt_x = red3::Grid1D::generate(uniform_x);
  REQUIRE(opt_x);
  red3::Grid1D x = opt_x.value();
  REQUIRE(5 == x.size());
  REQUIRE(0.25 == x.delta(0));
  // Check grid
  REQUIRE(0.0 == x[0]);
  REQUIRE(0.25 == x[1]);
  REQUIRE(0.5 == x[2]);
  REQUIRE(0.75 == x[3]);
  REQUIRE(1.0 == x[4]);

  red3::Uniform1D uniform_y(4, 1.0);
  std::optional<red3::Grid1D> opt_y = red3::Grid1D::generate(uniform_y);
  REQUIRE(opt_y);
  red3::Grid1D y = opt_y.value();
  REQUIRE(5 == y.size());
  REQUIRE(0.25 == y.delta(0));
  // Check grid
  REQUIRE(0.0 == y[0]);
  REQUIRE(0.25 == y[1]);
  REQUIRE(0.5 == y[2]);
  REQUIRE(0.75 == y[3]);
  REQUIRE(1.0 == y[4]);
}

