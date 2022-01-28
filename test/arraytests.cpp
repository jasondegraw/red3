// Copyright (C) 2015-2022 Jason W. DeGraw
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
#include "ut-wrapper.hpp"

#include "red3.hpp"
#include "array.hpp"

struct ArrayParentU
{
  ArrayParentU(red3::index_t ni, red3::index_t nj, red3::index_t nk = 1) : nu(ni+1), nj(nj), nk(nk)
  { }

  red3::index_t nu;
  red3::index_t nj;
  red3::index_t nk;
};

struct ArrayParentV
{
  ArrayParentV(red3::index_t ni, red3::index_t nj, red3::index_t nk = 1) : ni(ni), nv(nj + 1), nk(nk)
  {
  }

  red3::index_t ni;
  red3::index_t nv;
  red3::index_t nk;
};

boost::ut::suite arrays = [] {
  using namespace boost::ut;

  "basic 2d U arrays"_test = [] {
    ArrayParentU grid(4, 5);
    red3::ChildArrayU<ArrayParentU> array(&grid);
    expect(eq(25, array.size()) >> fatal);
    for (auto j = 0; j < grid.nj; j++) {
      for (auto i = 0; i < grid.nu; i++) {
        array(i, j, 0) = i;
      }
    }
    std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0 };
    for (auto i = 0; i < 25; i++) {
      expect(eq(v[i], array[i])) << "Index: " << i;
      array[i] = (double)i;
    }
    double value = 0.0;
    for (auto i = 0; i < 25; i++) {
      expect(eq(value, array[i])) << "Index: " << i;
      value += 1.0;
    }

  };

  "ops on 2d U arrays"_test = [] {
    ArrayParentU grid(4, 5);
    red3::ChildArrayU<ArrayParentU> array0(&grid);
    expect(eq(25, array0.size()) >> fatal);
    for (auto j = 0; j < grid.nj; j++) {
      for (auto i = 0; i < grid.nu; i++) {
        array0(i, j, 0) = i;
      }
    }

    red3::ChildArrayU<ArrayParentU> array1;
    expect(eq(array1.parent(), nullptr)); // Maybe this isn't the greatest idea

    array1 = array0;
    expect(eq(array0.parent(), array1.parent())); // Maybe this isn't the greatest idea either

    std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 2.0, 3.0, 4.0 };
    for (auto i = 0; i < 25; i++) {
      expect(eq(v[i], array1[i])) << "Index: " << i;
      array1[i] = (double)i;
    }
    double value = 0.0;
    for (auto i = 0; i < 25; i++) {
      expect(eq(value, array1[i])) << "Index: " << i;
      value += 1.0;
    }

    red3::ChildArrayU<ArrayParentU> array2;
    expect(eq(array2.parent(), nullptr));
    expect(array1 != array2) << "array1 != array2";
    expect(array1 == array0) << "array1 == array0";

    red3::ChildArrayU<ArrayParentU> array3(&grid);
    expect(eq(array3.parent(), &grid));
    expect(eq(25, array3.size()) >> fatal);
    expect(array3 != array0) << "array3 != array0";
    expect(array3 != array1) << "array3 != array1";
    expect(array3 != array2) << "array3 != array2";

    for (auto i = 0; i < 25; i++) {
      expect(0.0_d == array3[i]) << "Index: " << i;
    }

    array3.swap(array0);
    expect(array3 != array0) << "array3 != array0";
    expect(array3 == array1) << "array3 == array1";
    expect(array3 != array2) << "array3 != array2";

    value = 0.0;
    for (auto i = 0; i < 25; i++) {
      expect(0.0_d == array0[i]) << "Index: " << i;
      expect(eq(value, array1[i])) << "Index: " << i;
      expect(eq(value, array3[i])) << "Index: " << i;
      value += 1.0;
    }
  };

  "basic 2d V arrays"_test = [] {
    ArrayParentV grid(4, 5);
    red3::ChildArrayV<ArrayParentV> array(&grid);
    expect(eq(24, array.size()) >> fatal);
    for (auto j = 0; j < grid.nv; j++) {
      for (auto i = 0; i < grid.ni; i++) {
        array(i, j, 0) = j;
      }
    }
    std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    for (auto i = 0; i < 24; i++) {
      expect(eq(v[i], array[i])) << "Index: " << i;
      array[i] = (double)i;
    }
    double value = 0.0;
    for (auto i = 0; i < 24; i++) {
      expect(eq(value, array[i])) << "Index: " << i;
      value += 1.0;
    }

  };

  "ops on 2d V arrays"_test = [] {
    ArrayParentV grid(4, 5);
    red3::ChildArrayV<ArrayParentV> array0(&grid);
    expect(eq(24, array0.size()) >> fatal);
    for (auto j = 0; j < grid.nv; j++) {
      for (auto i = 0; i < grid.ni; i++) {
        array0(i, j, 0) = j;
      }
    }

    red3::ChildArrayV<ArrayParentV> array1;
    expect(eq(array1.parent(), nullptr)); // Maybe this isn't the greatest idea

    array1 = array0;
    expect(eq(array0.parent(), array1.parent())); // Maybe this isn't the greatest idea either

    std::vector<double> v = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0,
      5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    for (auto i = 0; i < 24; i++) {
      expect(eq(v[i], array1[i])) << "Index: " << i;
      array1[i] = (double)i;
    }
    double value = 0.0;
    for (auto i = 0; i < 24; i++) {
      expect(eq(value, array1[i])) << "Index: " << i;
      value += 1.0;
    }

    red3::ChildArrayV<ArrayParentV> array2;
    expect(eq(array2.parent(), nullptr));
    expect(array1 != array2) << "array1 != array2";
    expect(array1 == array0) << "array1 == array0";

    red3::ChildArrayV<ArrayParentV> array3(&grid);
    expect(eq(array3.parent(), &grid));
    expect(eq(24, array3.size()) >> fatal);
    expect(array3 != array0) << "array3 != array0";
    expect(array3 != array1) << "array3 != array1";
    expect(array3 != array2) << "array3 != array2";

    for (auto i = 0; i < 24; i++) {
      expect(0.0_d == array3[i]) << "Index: " << i;
    }

    array3.swap(array0);
    expect(array3 != array0) << "array3 != array0";
    expect(array3 == array1) << "array3 == array1";
    expect(array3 != array2) << "array3 != array2";

    value = 0.0;
    for (auto i = 0; i < 24; i++) {
      expect(0.0_d == array0[i]) << "Index: " << i;
      expect(eq(value, array1[i])) << "Index: " << i;
      expect(eq(value, array3[i])) << "Index: " << i;
      value += 1.0;
    }
  };

};


/*
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
*/

