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
#ifdef __cpp_modules
#ifdef _MSC_VER
#define BOOST_UT_DISABLE_MODULE
#include <boost/ut.hpp>
#else
import boost.ut; // Doesn't appear to work yet with MSVC/CMake
#endif
#else
#include <boost/ut.hpp>
#endif

#ifdef _MSC_VER
#include <windows.h>
#endif

int main(int argc, const char** argv)
{
  using namespace boost::ut;

  bool dry_run{false};

#ifdef _MSC_VER
  auto custom_colors = colors{ .none = "", .pass = "", .fail = "" };
  auto console_handle = GetStdHandle(STD_OUTPUT_HANDLE);
  if (console_handle != INVALID_HANDLE_VALUE) {
    DWORD console_mode;
    GetConsoleMode(console_handle, &console_mode);
    console_mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    auto success = SetConsoleMode(console_handle, console_mode);
    if (success != 0) {
      custom_colors.none = "\033[0m";
      custom_colors.pass = "\033[32m";
      custom_colors.fail = "\033[91m";
    }
  }
#else
  auto custom_colors = colors{};
#endif

  cfg<override> = {.filter = argc > 1 ? argv[1] : "",
                   .colors = custom_colors,
                   .dry_run = dry_run};

  return EXIT_SUCCESS;
}