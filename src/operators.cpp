#include "red3.h"

void divg(real_t *g, real_t *u, real_t *v, real_t *w, real_t *dx, 
          real_t *dy, real_t *dz, int nx, int ny, int nz)
{
  int i,j,k;
  int nxp1 = nx+1;
  int di = nx+1;
  int dj = 1;
  int ydelta = nx+2;
  int ip,im;
  int jp,jm;
  int ix,iy,iz;
  int n = 0;
  if(nz==1)
    {
      for(j=0;j<ny;j++)
        {
          im = n+di;
          jm = n+dj;
          for(i=0;i<nx;i++)
            {
              g[n] = (u[im+1]-u[im])/dx[i]
                +    (v[jm+ydelta]-v[jm])/dy[i];
              im++;
              jm++;
            }
          dj += 2;
          di += 1;
        }
    }
  else
    {
      // 3D goes here
    }
}
