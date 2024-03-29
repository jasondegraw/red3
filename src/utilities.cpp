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
#include "red3/utilities.hpp"
#include <stdexcept>
#include <stdlib.h>

namespace red3 {

double power2(double v)
{
  return v*v;
}

template <> double power<2>(double v)
{
  return v*v;
}

template <> double power<1>(double v)
{
  return v;
}

template <> double power<0>(double v)
{
  return 1.0;
}

void fatal(const std::string& mesg)
{
  throw std::runtime_error("red3: " + mesg);
}

void* callocate(size_t num, size_t size, const std::string& name)
{
  void* v = calloc(num, size);
  if (v == 0) {
    fatal("Failed to allocate " + name);
  }
  return v;
}

}

