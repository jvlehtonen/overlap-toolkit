/*
 * Copyright (c) 2022-2023 Jukka V. Lehtonen
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

#include <iostream>
#include <iomanip>

#include "Point.h"

Point& Point::operator+= ( const Point & rhs )
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

Point operator+ ( Point lhs, const Point & rhs )
{
  return lhs += rhs;
}

Point& Point::operator-= ( const Point & rhs )
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}

Point operator- ( Point lhs, const Point & rhs )
{
  return lhs -= rhs;
}

Point operator/ ( Point lhs, double rhs )
{
  lhs.x /= rhs;
  lhs.y /= rhs;
  lhs.z /= rhs;
  return lhs;
}

double dot( const Point & lhs, const Point & rhs )
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

std::ostream& operator<< ( std::ostream & out, const Point & pos )
{
  out << std::setw(9) << pos.x << ' ' << std::setw(9) << pos.y << ' ' << std::setw(9) << pos.z;
  return out;
}
