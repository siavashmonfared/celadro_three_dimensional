/*
 * This file is part of CELADRO_3D, Copyright (C) 2019-2021, Siavash Monfared
 *
 * This program is free software: you can redistribute it and/or modify
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
 
#ifndef STENCIL_HPP_
#define STENCIL_HPP_

/** Indices of a given point and its direct neighbourhood
 *
 * A stencil is a list of 9 (in two dimensions) grid indices of a point and its
 * neighbourhood. The convention is the following: if s is a given stencil, then
 * s[0][0] is the grid coordinate of the central point, s[1][0] the coordinate
 * of the point to its right, s[-1][0] the point to its left, etc.
 * */


struct stencil
{
  struct shifted_array_a
  {
    unsigned data[3];

    // CUDA_host_device
    unsigned& operator[](int i)
    { return data[i+1]; }

    // CUDA_host_device
    const unsigned& operator[](int i) const
    { return data[i+1]; }
  };

  struct shifted_array_b
  {
  shifted_array_a data[3];

  // CUDA_host_device
  shifted_array_a & operator[](int i)
  { return data[i+1]; }

  // CUDA_host_device
  const shifted_array_a & operator[](int i) const
  { return data[i+1]; }
  };
  
  shifted_array_b data[3];

  // CUDA_host_device
  shifted_array_b& operator[](int i)
  { return data[i+1]; }

  // CUDA_host_device
  const shifted_array_b& operator[](int i) const
  { return data[i+1]; }
  
};





#endif//STENCIL_HPP_
