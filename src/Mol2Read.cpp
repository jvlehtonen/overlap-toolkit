/*
 * Copyright (c) 2017-2023 Jukka V. Lehtonen
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

/*****************************************************************************
 * Includes
 ****************************************************************************/
#include <vector>
#include <iostream>

#include <QtCore>

#include "Mol2Read.h"

std::ostream& operator<< ( std::ostream& out, const Molecule& mol )
{
    out << "@<TRIPOS>MOLECULE\n";
    out << qPrintable(mol.name) << '\n';

    out << mol.atoms.size();
    if ( 0 < mol.bonds.size() )
    {
      out << ' ' << mol.bonds.size();

      if ( 0 < mol.substructures.size() ) out << ' ' << mol.substructures.size();
    }
    out << '\n';
    out << "SMALL\n";
    out << "USER_CHARGES\n";
    out << '\n';

    out << "@<TRIPOS>ATOM\n";
    unsigned long num {0};
    for ( auto atom : mol.atoms ) {
      atom.front() = QString::number( ++num );
      out << qPrintable( atom.join(' ') ) << '\n';
    }
    out << '\n';
    return out;
}


/****************************************************************************/
/*!
  \param Filename - name of file to parse.
*/
/****************************************************************************/
std::vector<Molecule>
parse( QTextStream & istr )
{
  std::vector<Molecule> result;
  QString              Mol_name;
  std::vector<QString> substructure;
  std::vector<QStringList> atoms;
  std::vector<QString> bonds;

  QString s;
  QStringList MolNumbers;
  unsigned long Line = 0;
  std::vector<QString>::size_type Atoms = 0;
  bool NextLine = true;

  // Read the file
  while ( !istr.atEnd() )
    {        // until end of file...
      if ( NextLine )
        {
          s = istr.readLine().simplified(); // line of text excluding '\n'
          auto index = s.indexOf( '#' );
          if ( -1 != index ) s = s.left( index ).simplified();
        }
      else
        {
          NextLine = true;
        }

      if ( ! s.isEmpty() )
        {
          if ( s.left( 17 ) == "@<TRIPOS>MOLECULE" )
            {
              // New molecule: create previous
              if ( 0 < atoms.size() )
                {
                  Q_ASSERT( Atoms == atoms.size() );
                  result.emplace_back( Mol_name, atoms, bonds, substructure );
                }

              // reset objects
              Mol_name.clear();
              substructure.clear();
              atoms.clear();
              bonds.clear();

              Atoms = 0;

              Line = 0;
              int CommentStarts = -1;
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  CommentStarts = s.indexOf( '#' );
                  if ( -1 != CommentStarts )
                    {
                      s = s.left( s.indexOf( '#' ) ).simplified();
                    }

                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  else if ( 0 != CommentStarts )
                    {
                      switch( Line )
                        {
                        case 0:
                          Mol_name = s;
                          if ( s.isEmpty() ) s = "Mol2";
                          qDebug() << "Mol_name: " << s << '\n';
                          ++Line;
                          break;
                        case 1:
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
                          MolNumbers = s.split( ' ', QString::SkipEmptyParts );
#else
                          MolNumbers = s.split( ' ', Qt::SkipEmptyParts );
#endif
                          Atoms = MolNumbers.first().toInt();
                          atoms.reserve( Atoms );
                          qDebug() << "Num_atoms: " << Atoms << '\n';
                          ++Line;
                          break;
                        case 2:
                          ++Line;
                          break;
                        case 3:
                          qDebug() << "Charge_type: " << s << '\n';
                          ++Line;
                          break;
                        default:
                          break;
                        }
                    }
                }
            }
          else if ( s.left( 13 ) == "@<TRIPOS>DICT" )
            {
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  if ( -1 != s.indexOf( '#' ) ) s = s.left( s.indexOf( '#' ) ).simplified();
                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  // We should do something here
                }
            }
          else if ( s.left( 13 ) == "@<TRIPOS>ATOM" )
            {
              Line = 0;
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  if ( -1 != s.indexOf( '#' ) ) s = s.left( s.indexOf( '#' ) ).simplified();
                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  else if ( ! s.isEmpty() )
                    {
                      if ( Line < Atoms )
                        {
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
                          auto entry = s.split( ' ', QString::SkipEmptyParts );
#else
                          auto entry = s.split( ' ', Qt::SkipEmptyParts );
#endif
                          atoms.push_back( entry );
                          ++Line;
                        }
                    }
                }
            }
          else if ( s.left( 13 ) == "@<TRIPOS>BOND" )
            {
              Line = 0;
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  if ( -1 != s.indexOf( '#' ) ) s = s.left( s.indexOf( '#' ) ).simplified();
                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  else if ( ! s.isEmpty() )
                    {
                      bonds.push_back( s );
                      ++Line;
                    }
                }
            }
          else if ( s.left( 21 ) == "@<TRIPOS>SUBSTRUCTURE" )
            {
              Line = 0;
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  if ( -1 != s.indexOf( '#' ) ) s = s.left( s.indexOf( '#' ) ).simplified();
                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  else if ( ! s.isEmpty() )
                    {
                      substructure.push_back( s );
                      ++Line;
                    }
                }
            }
          else if ( s.left( 12 ) == "@<TRIPOS>SET" )
            {
              while ( !istr.atEnd() )
                {
                  s = istr.readLine().simplified();
                  if ( -1 != s.indexOf( '#' ) ) s = s.left( s.indexOf( '#' ) ).simplified();
                  if ( s.left( 9 ) == "@<TRIPOS>" )
                    {
                      NextLine = false;
                      break;
                    }
                  // We should do something here
                }
            }
        }
    }

  // Create last molecule
  Q_ASSERT( Atoms == atoms.size() );
  if ( ! atoms.empty() )
    {
      result.emplace_back( Mol_name, atoms, bonds, substructure );
    }

  return result;
}
