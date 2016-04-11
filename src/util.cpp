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
#include "util.hpp"
#include <stdexcept>
#include <stdlib.h>

namespace red3 {

void fatal(const std::string &mesg)
{
  throw std::runtime_error("red3: " + mesg);
}

void* callocate(size_t num, size_t size, const std::string &name)
{
  void* v = calloc(num, size);
  if(v == 0) {
    fatal("Failed to allocate " + name);
  }
  return v;
}

}
