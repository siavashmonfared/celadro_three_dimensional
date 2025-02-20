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


#ifndef FILES_HPP_
#define FILES_HPP_

/** Create directory
 *
 * All subdirectories are created as necessary.
 * */
void create_directory(const std::string& dir);

/** Remove file or directory
 *
 * This is equivalent to 'rm -rf fname'
 * */
void remove_file(const std::string& fname);

/** Compress file iname to oname.zip using zip
 *
 * The file is moved to the archive, i.e. iname does not exsit after this
 * function has returned.
 * */
void compress_file(const std::string& iname, const std::string& oname);

#endif//FILES_HPP_
