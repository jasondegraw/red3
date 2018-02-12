// Copyright (C) 2015-2017 Jason W. DeGraw
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
#include "upwind.hpp"
#include <stdexcept>
#include <stdlib.h>
#include "arrayops.hpp"

namespace red3 {
namespace upwind {

void StaggeredIncompressibleSteadyFlow::setupU()
{
  // Main body of block
  //
  //              N
  //              |
  //              |
  // +-----N1-----+----N2----+
  // |      ::::::|:::::     |
  // W     p1:::::P::::p2    E
  // |      ::::::|:::::     |
  // +-----S1-----+----S2----+
  //              |
  //              S
  //
  int k = 0;
  for (int j = 1; j < nj - 1; j++) {
    for (int i = 1; i < ni; i++) {
      int ijk(UINDEX(i, j, k, nu, nj, nk));
      // East
      int east(ijk + 1);
      // West
      int west(ijk - 1);
      // South
      int south2(VINDEX(i, j, k, ni, nv, nk));
      int north2(south2 + 1);
      // North
      int south1(south2 - nv);
      int north1(south1 + 1);

      // Compute fluxes
      double Fe = rho*0.5*(u(east) + u(ijk))*dy[j];
      double Fw = rho*0.5*(u(ijk) + u(west))*dy[j];
      double Fn = rho*0.5*(v(north1)*dx[i - 1] + v(north2)*dx[i]);
      double Fs = rho*0.5*(v(south1)*dx[i - 1] + v(south2)*dx[i]);

      double De = mu*dy[j] / dx[i];
      double Dw = mu*dy[j] / dx[i - 1];
      double Dn = mu*(dx[i - 1] + dx[i]) / (dy[j] + dy[j + 1]); // 0.5s cancel out
      double Ds = mu*(dx[i - 1] + dx[i]) / (dy[j - 1] + dy[j]); // 0.5s cancel out

      // Compute the coefficients
      ae[ijk] = De*A(std::abs(Fe / De)) + std::max(-Fe, 0.0);
      aw[ijk] = Dw*A(std::abs(Fw / Dw)) + std::max( Fw, 0.0);
      an[ijk] = Dn*A(std::abs(Fn / Dn)) + std::max(-Fn, 0.0);
      as[ijk] = Ds*A(std::abs(Fs / Ds)) + std::max( Fs, 0.0);

      double ap(ae[ijk] + aw[ijk] + an[ijk] + as[ijk]);

      // Pre-divide by the central coefficent
      ae[ijk] /= ap;
      aw[ijk] /= ap;
      an[ijk] /= ap;
      as[ijk] /= ap;

      // Compute the forcing term
      int p2(INDEX(i, j, k, ni, nj, nk));
      int p1(p2 - 1);
      b[ijk] = 2.0*(p(p1) - p(p2)) / (ap*(dx[i - 1] + dx[i]));
    }
  }
  // South boundary of block
  //
  //              N
  //              |
  //              |
  // +-----N1-----+----N2----+
  // |      ::::::|:::::     |
  // W     p1:::::P::::p2    E
  // |      ::::::|:::::     |
  // +-----S1-----b----S2----+
  //              |
  //              S
  //
  // Mirror cell method: require that the velocity at the boundary is the specified value U_b by computing
  //
  //   U_b = 0.5*(U_P + U_S)
  //
  // Thus we can calculate the mirror cell value U_S = 2*U_b - U_P
  k = 0;
  int j = 0;
  for (int i = 1; i < ni; i++) {
    int ijk(UINDEX(i, j, k, nu, nj, nk));
    // East
    int east(ijk + 1);
    // West
    int west(ijk - 1);
    // South
    int south2(VINDEX(i, j, k, ni, nv, nk));
    int north2(south2 + 1);
    // North
    int south1(south2 - nv);
    int north1(south1 + 1);

    // Compute fluxes
    double Fe = rho*0.5*(u(east) + u(ijk))*dy[j];
    double Fw = rho*0.5*(u(ijk) + u(west))*dy[j];
    double Fn = rho*0.5*(v(north1)*dx[i - 1] + v(north2)*dx[i]);
    double Fs = rho*0.5*(v(south1)*dx[i - 1] + v(south2)*dx[i]);

    double De = mu*dy[j] / dx[i];
    double Dw = mu*dy[j] / dx[i - 1];
    double Dn = mu*(dx[i - 1] + dx[i]) / (dy[j] + dy[j + 1]); // 0.5s cancel out
    // The mirror cell is the same size as cell j, so dy[j-1] = dy[j]
    // Ds = mu*(dx[i - 1] + dx[i]) / (dy[j - 1] + dy[j]) = mu*(dx[i - 1] + dx[i]) / (dy[j] + dy[j])
    double Ds = 0.5*mu*(dx[i - 1] + dx[i]) / dy[j]; // 2 0.5s cancel out, one remains

    // Compute the coefficients
    ae[ijk] = De*A(std::abs(Fe / De)) + std::max(-Fe, 0.0);
    aw[ijk] = Dw*A(std::abs(Fw / Dw)) + std::max( Fw, 0.0);
    an[ijk] = Dn*A(std::abs(Fn / Dn)) + std::max(-Fn, 0.0);
    as[ijk] = Ds*A(std::abs(Fs / Ds)) + std::max( Fs, 0.0);

    // Tack on an additional term from the "s" fluxes: as*U_S = as*(2*U_b - U_P) = 2*as*U_b - as*U_P
    // Since this added term is from the RHS, it's an added as
    double ap(ae[ijk] + aw[ijk] + an[ijk] + 2*as[ijk]);

    // Pre-divide by the central coefficent
    ae[ijk] /= ap;
    aw[ijk] /= ap;
    an[ijk] /= ap;
    as[ijk] /= ap;

    // Compute the forcing term
    int p2(INDEX(i, j, k, ni, nj, nk));
    int p1(p2 - 1);
    // Add the boundary u to the forcing term
    b[ijk] = 2*(p(p1) - p(p2)) / (ap*(dx[i - 1] + dx[i])) + as[ijk]*2*southU(i, k);

  }

  // North boundary of block
  //
  //              N
  //              |
  //              |
  // +-----N1-----+----N2----+
  // |      ::::::|:::::     |
  // W     p1:::::P::::p2    E
  // |      ::::::|:::::     |
  // +-----S1-----b----S2----+
  //              |
  //              S
  //
  // Mirror cell method: require that the velocity at the boundary is the specified value U_b by computing
  //
  //   U_b = 0.5*(U_P + U_N)
  //
  // Thus we can calculate the mirror cell value U_N = 2*U_b - U_P
  k = 0;
  j = nj - 1;
  for (int i = 1; i < ni; i++) {
    int ijk(UINDEX(i, j, k, nu, nj, nk));
    // East
    int east(ijk + 1);
    // West
    int west(ijk - 1);
    // South
    int south2(VINDEX(i, j, k, ni, nv, nk));
    int north2(south2 + 1);
    // North
    int south1(south2 - nv);
    int north1(south1 + 1);

    // Compute fluxes
    double Fe = rho*0.5*(u(east) + u(ijk))*dy[j];
    double Fw = rho*0.5*(u(ijk) + u(west))*dy[j];
    double Fn = rho*0.5*(v(north1)*dx[i - 1] + v(north2)*dx[i]);
    double Fs = rho*0.5*(v(south1)*dx[i - 1] + v(south2)*dx[i]);

    double De = mu*dy[j] / dx[i];
    double Dw = mu*dy[j] / dx[i - 1];
    // The mirror cell is the same size as cell j, so dy[j + 1] = dy[j]
    // Dn = mu*(dx[i - 1] + dx[i]) / (dy[j] + dy[j + 1]) = mu*(dx[i - 1] + dx[i]) / (dy[j] + dy[j])
    double Dn = 0.5*mu*(dx[i - 1] + dx[i]) / dy[j]; // 2 0.5s cancel out, one remains
    double Ds = mu*(dx[i - 1] + dx[i]) / (dy[j - 1] + dy[j]);  // 0.5s cancel out

    // Compute the coefficients
    ae[ijk] = De*A(std::abs(Fe / De)) + std::max(-Fe, 0.0);
    aw[ijk] = Dw*A(std::abs(Fw / Dw)) + std::max( Fw, 0.0);
    an[ijk] = Dn*A(std::abs(Fn / Dn)) + std::max(-Fn, 0.0);
    as[ijk] = Ds*A(std::abs(Fs / Ds)) + std::max( Fs, 0.0);

    // Tack on an additional term from the "n" fluxes: as*U_N = an*(2*U_b - U_P) = 2*an*U_b - an*U_P
    // Since this added term is from the RHS, it's an added an
    double ap(ae[ijk] + aw[ijk] + 2*an[ijk] + as[ijk]);

    // Pre-divide by the central coefficent
    ae[ijk] /= ap;
    aw[ijk] /= ap;
    an[ijk] /= ap;
    as[ijk] /= ap;

    // Compute the forcing term
    int p2(INDEX(i, j, k, ni, nj, nk));
    int p1(p2 - 1);
    // Add the boundary u to the forcing term
    b[ijk] = 2*(p(p1) - p(p2)) / (ap*(dx[i - 1] + dx[i])) + as[ijk] * 2*northU(i, k);

  }

}

}

}
