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

#include "staggeredgrid.hpp"
#include "gtest/gtest.h"
#include <iostream>

TEST(StaggerdGridTests, Basic2D)
{
  red3::StaggeredGrid grid(4, 5);
  EXPECT_EQ(4, grid.ni);
  EXPECT_EQ(5, grid.nu);
  EXPECT_EQ(5, grid.nj);
  EXPECT_EQ(6, grid.nv);
  EXPECT_EQ(1, grid.nk);
  EXPECT_EQ(0, grid.nw);
  // Check grid
  ASSERT_EQ(5, grid.x.size());
  std::vector<double> x = { 0.0, 0.25, 0.5, 0.75, 1.0 };
  for(unsigned i = 0; i < grid.nu; i++) {
    EXPECT_DOUBLE_EQ(x[i], grid.x[i]);
  }
  ASSERT_EQ(4, grid.xm.size());
  std::vector<double> xm = { 0.125, 0.375, 0.625, 0.875 };
  for(unsigned i = 0; i < grid.ni; i++) {
    EXPECT_DOUBLE_EQ(xm[i], grid.xm[i]);
  }
  ASSERT_EQ(6, grid.y.size());
  std::vector<double> y = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
  for(unsigned i = 0; i < grid.nv; i++) {
    EXPECT_DOUBLE_EQ(y[i], grid.y[i]);
  }
  ASSERT_EQ(5, grid.ym.size());
  std::vector<double> ym = { 0.1, 0.3, 0.5, 0.7, 0.9 };
  for(unsigned i = 0; i < grid.ni; i++) {
    EXPECT_DOUBLE_EQ(ym[i], grid.ym[i]);
  }
  // Check velocity
  grid.setU([](double x, double y){return x;});
  for(unsigned j = 0; j < grid.nj; j++) {
    for(unsigned i = 0; i < grid.nu; i++) {
      EXPECT_EQ(x[i], grid.u(i, j, 0));
    }
  }
}

