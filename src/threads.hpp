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

#ifndef THREADS_HPP_
#define THREADS_HPP_

#ifdef _OPENMP
  #define PRAGMA_OMP(x) _Pragma(#x)
#else
  #define PRAGMA_OMP(cmd) {}
#endif

void SetThreads();

#endif//THREADS_HPP_
