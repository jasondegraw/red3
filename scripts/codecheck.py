# Copyright (C) 2015-2016 Jason W. DeGraw
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Check code files for required elements

@author: Jason W. DeGraw
'''
import glob

notice0 = """# Copyright (C) 2015-2016 Jason W. DeGraw
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

notice1 = """// Copyright (C) 2015-2016 Jason W. DeGraw
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
"""

def check(files, text):
    count = 0
    has_license = 0
    missing_license = 0
    for filename in files:
        count += 1
        fp = open(filename,'r')
        txt = fp.read()
        fp.close()
        if not txt.startswith(text):
            print('%s does not have the proper license' % filename)
            missing_license += 1
        else:
            has_license += 1
    print('Checked %d files' % count)
    print('\t%d files passed' % has_license)
    print('\t%d files failed' % missing_license)
    
files = glob.iglob('../**/*.py')
check(files, notice0)

files = []
dirs = ['src', 'test', 'laminar-channel']
for directory in dirs:
    files.extend(glob.iglob('../%s/*.cpp'%directory))
    files.extend(glob.iglob('../%s/*.hpp'%directory))
check(files, notice1)

