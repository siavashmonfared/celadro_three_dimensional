/*
 * This file is part of CELADRO-3D, (C) 2020-2021, Siavash Monfared
 * and Romain Mueller (CELADRO). This program is free software: you can 
 * redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "header.hpp"
#include "model.hpp"
#include "derivatives.hpp"

using namespace std;

void Model::AddCellAtNode(unsigned n, unsigned q, const coord& center)
{
    
  const auto      k = GetIndexFromPatch(n, q);
  const unsigned xk = GetXPosition(k);
  const unsigned yk = GetYPosition(k);
  const unsigned zk = GetZPosition(k);

  const auto radius = max(R/2., 4.);
  // creates spheres 
  
  if(
      (BC==0 and pow(wrap(diff(yk, center[1]), Size[1]), 2)
       + pow(wrap(diff(xk, center[0]), Size[0]), 2) + pow(wrap(diff(zk, center[2]), Size[2]), 2) <=ceil(radius*radius))
      or
      (BC>=1 and pow(diff(yk, center[1]), 2)
       + pow(diff(xk, center[0]), 2) + pow(diff(zk, center[2]), 2) <=ceil(radius*radius))
    )

    
  {
            
    phi[n][q]     = 1.;
    phi_old[n][q] = 1.;
    vol[n]      += 1.;
    square[k]    += 1.;
    sum[k]       += 1.;
    sumQ00[k]    += Q00[n];
    sumQ01[k]    += Q01[n];    
  }
  else
  {
    phi[n][q]     = 0.;
    phi_old[n][q] = 0.;
  }
  
}

void Model::AddCell(unsigned n, const coord& center)
{
  // update patch coordinates
  patch_min[n] = (center+Size-patch_margin)%Size;
  patch_max[n] = (center+patch_margin-1u)%Size;

  // init polarisation and nematic
  theta_pol[n] = noise*Pi*(1-2*random_real());
  polarization[n] = { Spol*cos(theta_pol[n]), Spol*sin(theta_pol[n]) };
  
  // need to be developed for 3D (siavash 03/31/2021)
  theta_nem[n] = noise*Pi*(1-2*random_real());
  Q00[n] = Snem*cos(2*theta_nem[n]);
  Q01[n] = Snem*sin(2*theta_nem[n]);

  // create the cells at the centers we just computed
  for(unsigned q=0; q<patch_N; ++q)
    AddCellAtNode(n, q, center);
  
  com[n]   = vec<double, 3>(center);

}

void Model::AddCellMix(unsigned n, const coord& center)
{
  // update patch coordinates
  patch_min[n] = (center+Size-patch_margin)%Size;
  patch_max[n] = (center+patch_margin-1u)%Size;

  // init polarisation and nematic
  theta_pol[n] = noise*Pi*(1-2*random_real());
  polarization[n] = { Spol*cos(theta_pol[n]), Spol*sin(theta_pol[n]) };
  theta_nem[n] = noise*Pi*(1-2*random_real());
  Q00[n] = Snem*cos(2*theta_nem[n]);
  Q01[n] = Snem*sin(2*theta_nem[n]);

  // create the cells at the centers we just computed
  for(unsigned q=0; q<patch_N; ++q)
    AddCellAtNode(n, q, center);

  com[n]   = vec<double, 3>(center);
  
}

void Model::Configure()
{

  
  if(init_config=="input const")
  {
    string fname = "input_str.dat";
    double xcoor, ycoor, zcoor;
    fstream file(fname);
    for (unsigned n = 0 ; n < nphases ; n++){
    file >> xcoor >> ycoor >> zcoor;
    AddCellMix(n,{xcoor,ycoor,zcoor});// reads file for positions only 
    }
  }  


  else throw error_msg("error: initial configuration '",
      init_config, "' unknown.");
      
}

void Model::ConfigureWalls(int BC_)
{
  switch(BC_) // note that BC_ is argument, and different from class member BC
  {
  case 0:
    // no walls (pbc)
    for(unsigned k=0; k<N; ++k)
      walls[k] = 0;
    break;
      
  case 1:
    // Exponentially falling phase-field:
    for(unsigned k=0; k<N; ++k)
    {
      const double x = GetXPosition(k);
      const double y = GetYPosition(k);
      const double z = GetZPosition(k);
      // this is the easiest way: each wall contributes as an exponentially
      // falling potential and we do not care about overalps
      walls[k] =   exp(-y/wall_thickness)
                 + exp(-x/wall_thickness)
                 + exp(-z/wall_thickness)
                 + exp(-(Size[0]-1-x)/wall_thickness)
                 + exp(-(Size[1]-1-y)/wall_thickness)
                 + exp(-(Size[2]-1-z)/wall_thickness);
    }
    break;
  
  case 2:
    // constant wall
    for(unsigned k=0; k<N; ++k)
    {
      const double z = GetZPosition(k);
      
      if (z < wall_thickness) walls[k] = 1.;
      else walls[k] = 0.;
      
    }
    break;
      
  case 3:
    // box
    for(unsigned k=0; k<N; ++k)
    {
      const double x = GetXPosition(k);
      const double y = GetYPosition(k);
      const double z = GetZPosition(k);
      
      if (z < wall_thickness)             	walls[k] = 1.;
      else if(z > Size[2]-wall_thickness) 	walls[k] = 1.;
      else if(x < wall_thickness)         	walls[k] = 1.;
      else if(x > Size[0]-wall_thickness) 	walls[k] = 1.;
      else if(y < wall_thickness)         	walls[k] = 1.;
      else if(y > Size[1]-wall_thickness) 	walls[k] = 1.;
      else walls[k] = 0.;
      
    }
    break;
    
  
    default:
    throw error_msg("boundary condition unknown.");
  }
  // pre-compute derivatives
  for(unsigned k=0; k<N; ++k)
  {
    const auto& s = neighbors[k];

    walls_dx[k] = derivX(walls, s);
    walls_dy[k] = derivY(walls, s);
    walls_dz[k] = derivZ(walls, s);
    walls_laplace[k] = laplacian(walls, s);
  }
}


















