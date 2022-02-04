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
#ifndef RED3_HPP
#define RED3_HPP

#if _WIN32 || _MSC_VER

#ifdef red3_core_EXPORTS
#define RED3_API __declspec(dllexport)
#else
#define RED3_API __declspec(dllimport)
#endif
#else
#define RED3_API
#endif

#ifdef _MSC_VER
#pragma warning(disable: 4251)
#endif

namespace red3 {
using index_t = size_t;
}


#endif
