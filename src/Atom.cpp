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

#include <map>
#include <numeric>
#include <iostream>
#include <iomanip>

#include <QStringList>

#include "Atom.h"

Point Atom::pos() const
{
  auto avg = std::accumulate( begin(posit), end(posit), Point() );
  return posit.size() ? avg / posit.size() : Point() ;
}

std::ostream& print( std::ostream & out, const Atom & atom, unsigned long num )
{
  static std::map<QString,int> counters;
  auto type = atom.type.split(".").front();
  ++(counters[ type ]);
  out << std::fixed;
  out << std::setw(7) << num << ' ';
  out << std::left << std::setw(7) << qPrintable( QString("%1%2").arg(type).arg(counters[type]) );
  out.precision(4);
  out << std::right << ' ' << atom.pos();
  out << std::left << ' ' << std::setw(9) << qPrintable( atom.type );
  out.precision(3);
  out << std::right << " 1 LIG     " << std::setw(9) << atom.charge;
  return out;
}
