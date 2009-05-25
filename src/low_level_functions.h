/* Copyright (C) 2009 Carlo de Falco

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

int findspan(int n, int p, double u, const RowVector& U);
void basisfun(int i, double u, int p, const RowVector& U, RowVector& N);
void basisfunder (int i, int pl, double uu, const RowVector& u_knotl, int nders, NDArray& dersv);
