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

#ifndef atom_h
#define atom_h

#include <vector>
#include <iosfwd>
#include <QString>

#include "Point.h"

struct Atom {
  QString serial;
  QString name;
  std::vector<Point> posit;
  QString type;
  double  charge {};
  bool    mark   {false};

  Atom( const QString& serial, const QString& name, const Point& pos, const QString& type, double charge )
    : serial{serial}, name{name}, posit{pos}, type{type}, charge{charge}
  { }

  Point pos() const;
  bool mono() const { return 1 == posit.size(); }
};

std::ostream& print( std::ostream & out, const Atom & atom, unsigned long num );

#endif
