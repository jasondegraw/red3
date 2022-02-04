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

#define INDEX(i,j,k,ni,nj,nk) ((i)+(j)*(ni)+(k)*(ni)*(nj))
#define UINDEX(i,j,k,nu,nj,nk) ((i)+(j)*(nu)+(k)*(nu)*(nj))
#define VINDEX(i,j,k,ni,nv,nk) ((i)*(nv)+(j)+(k)*(ni)*(nv))
#define WINDEX(i,j,k,ni,nj,nw) ((i)*(nw)+(j)*(ni)*(nw)+(k))
#define IKINDEX(i,k,ni,nk) ((i)+(k)*(ni))
