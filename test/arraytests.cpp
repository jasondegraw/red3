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

#include "red3.hpp"
#include "array.hpp"
#include "catch.hpp"

struct ArrayParentU
{
  ArrayParentU(red3::index_t ni, red3::index_t nj, red3::index_t nk = 1) : nu(ni+1), nj(nj), nk(nk)
  { }

  red3::index_t nu;
  red3::index_t nj;
  red3::index_t nk;
};

TEST_CASE("Basic 2D U Arrays", "[array]")
{
  ArrayParentU grid(4, 5);
  red3::ChildArrayU<ArrayParentU> array(&grid);
  REQUIRE(25 == array.size());
  for(auto j = 0; j < grid.nj; j++) {
    for(auto i = 0; i < grid.nu; i++) {
      array(i, j, 0) = i;
    }
  }
  std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0,
    0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0 };
  for(auto i = 0; i < 25; i++) {
    INFO("Index: " << i);
    REQUIRE(v[i] == array[i]);
    array[i] = (double)i;
  }
  double value = 0.0;
  for(auto i = 0; i < 25; i++) {
    INFO("Index: " << i);
    REQUIRE(value == array[i]);
    value += 1.0;
  }
}

struct ArrayParentV
{
  ArrayParentV(red3::index_t ni, red3::index_t nj, red3::index_t nk = 1) : ni(ni), nv(nj+1), nk(nk)
  {
  }

  red3::index_t ni;
  red3::index_t nv;
  red3::index_t nk;
};

TEST_CASE("Basic 2D V Arrays", "[array]")
{
  ArrayParentV grid(4, 5);
  red3::ChildArrayV<ArrayParentV> array(&grid);
  REQUIRE(24 == array.size());
  for (auto j = 0; j < grid.nv; j++) {
    for (auto i = 0; i < grid.ni; i++) {
      array(i, j, 0) = j;
    }
  }
  std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0,
    5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  for (auto i = 0; i < 24; i++) {
    INFO("Index: " << i);
    REQUIRE(v[i] == array[i]);
    array[i] = (double)i;
  }
  double value = 0.0;
  for (auto i = 0; i < 24; i++) {
    INFO("Index: " << i);
    REQUIRE(value == array[i]);
    value += 1.0;
  }
}
