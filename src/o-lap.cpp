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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

#include <QtCore>

#include "Mol2Read.h"
#include "Point.h"
#include "Atom.h"
#include "json.h"

void header( std::ostream& out, const QString& name, size_t atoms )
{
    out << "@<TRIPOS>MOLECULE\n";
    out << ' ' << qPrintable(name) << '\n';
    out << ' ' << atoms << '\n';
    out << " SMALL\n";
    out << " USER_CHARGES\n";
    out << '\n';

    out << "@<TRIPOS>ATOM\n";
}

void merge( std::ostream& out, const std::vector<Molecule>& mols, const std::vector<bool>& skipped )
{
    size_t atoms {0};
    for ( size_t i=0; i < mols.size(); ++i ) {
      if ( not skipped[i] ) {
        atoms += mols[i].atoms.size();
      }
    }
    header( out, "Merged fragments", atoms );
    size_t num {0};
    for ( size_t i=0; i < mols.size(); ++i ) {
      if ( not skipped[i] ) {
        for ( auto atom : mols[i].atoms ) {
          atom.front() = QString::number( ++num );
          out << qPrintable( atom.join(' ') ) << '\n';
        }
      }
    }
    out << '\n';
}

std::map<QString,int> atomtypes
{
};


//!
//! Show atomtypes as categories and in JSON
//!
void showsimilar()
{
  QMap<QString, QVariant> smap;
  std::map<int,std::vector<QString>> typecats;
  for ( const auto& type : atomtypes ) {
    typecats[type.second].emplace_back( type.first );
    smap.insert( type.first, type.second );
  }
  std::cout << std::string( 30, '#' ) << '\n';
  std::cout << "# Category: types\n";
  std::cout << std::string( 30, '#' ) << '\n';
  for ( const auto& cat : typecats ) {
    std::cout << cat.first << ':';
    for ( const auto& type : cat.second ) {
      std::cout << ' ' << qPrintable( type );
    }
    std::cout << '\n';
  }
  std::cout << std::string( 30, '#' ) << '\n';
  std::cout << "# in JSON:\n";
  std::cout << std::string( 30, '#' ) << '\n';
  auto obj = QJsonObject::fromVariantMap( smap );
  QJsonDocument doc { obj };
  auto arr = doc.toJson( QJsonDocument::Indented );
  std::cout << qPrintable( arr );
}


int atomtype( const QString & lhs )
{
  auto lp = atomtypes.find( lhs );
  if ( lp != atomtypes.end() ) return lp->second;
  return 0;
}

bool sametype( const Atom & lhs, const Atom & rhs, bool similar, double charge = 0.2 )
{
  if ( charge < std::abs(lhs.charge - rhs.charge) ) return false;
  if ( similar ) {
    auto lp = atomtypes.find( lhs.type );
    auto rp = atomtypes.find( rhs.type );
    if ( lp != atomtypes.end() && rp != atomtypes.end() ) return lp->second == rp->second;
  }
  return lhs.type == rhs.type;
}

double sqrdist( const Atom & lhs, const Atom & rhs )
{
  auto pos = lhs.pos() - rhs.pos();
  auto dist = sqrt( dot( pos, pos ) );
  if ( lhs.mono() && rhs.mono() ) dist *= 2;
  return (lhs.type == "C.ar" || lhs.type == "N.ar") && 1.38 < dist ? (4.0 * dist) : dist;
}


double distance( const Atom & lhs, const Atom & rhs )
{
  auto pos = lhs.pos() - rhs.pos();
  auto dist = sqrt( dot( pos, pos ) );
  if ( lhs.mono() && rhs.mono() ) dist *= 2;
  return dist;
}


double sdist( const Atom & lhs, const Atom & rhs )
{
  auto pos = lhs.pos() - rhs.pos();
  return dot( pos, pos );
}


//!
//! Fuse atoms based on MCL-style cluster data
//!
void mcl2atoms( QTextStream& istr, std::map<int,std::vector<Atom>>& atomcats,
                size_t molecule, QString prefix, int argc, char *argv[],
                unsigned long cmin, unsigned long cminchr, double nibthreshold )
{
//  const auto original_count = atoms.size();
  unsigned long count {0};
  QString line;
  while ( istr.readLineInto(&line) ) {
    ++count;
  }

  std::set<std::pair<unsigned long,unsigned long>> used;
  if ( 0 < count && istr.seek( 0 ) ) {
    unsigned long serial {0};
    std::ostringstream ostr;
    while ( istr.readLineInto(&line) ) {
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
      auto words = line.split('\t', QString::SkipEmptyParts);
#else
      auto words = line.split("\t", Qt::SkipEmptyParts);
#endif
      if ( 0 < words.size() ) {
        auto id = words[0].split( "_" );
        unsigned long tnum = id[0].toULong();
        unsigned long anum = id[1].toULong();
        used.insert( std::make_pair( tnum, anum ) );
        auto& to = atomcats.at( tnum )[ anum ].posit;
        double& toc = atomcats.at( tnum )[ anum ].charge;
        for ( int w = 1; w < words.size(); ++w ) {
          id = words[w].split( "_" );
          unsigned long wtnum = id[0].toULong();
          unsigned long wanum = id[1].toULong();
          used.insert( std::make_pair( wtnum, wanum ) );
          const auto& from = atomcats.at( wtnum )[ wanum ].posit;
          to.insert( to.end(), from.begin(), from.end() );

          const double frc = atomcats.at( wtnum )[ wanum ].charge;
          if ( std::abs(toc) < std::abs(frc) ) toc = frc;
        }
        if ( std::abs(toc) <= nibthreshold ) {
          if ( cmin <= to.size() ) {
            ++serial;
            print( ostr, atomcats.at( tnum )[ anum ], serial );
            ostr << '\n';
          }
        } else {
          if ( cminchr <= to.size() ) {
            ++serial;
            print( ostr, atomcats.at( tnum )[ anum ], serial );
            ostr << '\n';
          }
        }
      }
    }

    // if MCL output does not contain all single atom clusters
    // then must print the rest separately
    if ( cmin < 2 ) {
      for ( const auto& cat : atomcats ) {
        unsigned long tnum = cat.first;
        for ( unsigned long anum {}; anum < cat.second.size(); ++anum ) {
          if ( used.find( std::make_pair( tnum, anum )  ) == used.end() ) {
            if ( std::abs(cat.second[ anum ].charge) <= nibthreshold ) {
              if ( cmin < 2 ) {
                ++serial;
                print( ostr, cat.second[ anum ], serial );
                ostr << '\n';
              }
            } else {
              if ( cminchr <= 2 ) {
                ++serial;
                print( ostr, cat.second[ anum ], serial );
                ostr << '\n';
              }
            }
          }
        }
      }
    }

    std::cout << "# Output from overlap " << qPrintable( QCoreApplication::applicationVersion() ) << '\n';
    std::cout << "# Created: " << qPrintable(QDateTime::currentDateTime().toString()) << '\n';
    std::cout << "# Command:";
    for (int a{}; a < argc; ++a ) std::cout << ' ' << argv[a];
    std::cout << "\n\n";
    header( std::cout, QString("%1%2").arg(prefix).arg(molecule), serial );
    std::cout << ostr.str();
  }
}


//!
//! Show types of atoms in MCL-style cluster data
//!
void mcl2types( QTextStream& istr, std::map<int,std::vector<Atom>>& atomcats, QTextStream& ostr )
{
  QString line;
  while ( istr.readLineInto(&line) ) {
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    auto words = line.split('\t', QString::SkipEmptyParts);
#else
    auto words = line.split("\t", Qt::SkipEmptyParts);
#endif
    if ( 0 < words.size() ) {
      for ( auto w : words ) {
        const auto id = w.split( "_" );
        unsigned long tnum = id[0].toULong();
        unsigned long anum = id[1].toULong();
        ostr << QString("%1").arg( atomcats.at( tnum )[ anum ].type, -6 );
      }
      ostr << '\n';
    }
  }
}


void internal_merge( std::vector<Atom>& atoms, unsigned long cmin, unsigned long cminchr,
                     double nibthreshold, double cutoff,
                     bool similar, double chargediff, QMap<QString, QVariant> cutmap )
{
  double best = 99999999.0;
  while ( 1 < atoms.size() ) {
    std::vector<double> distmat( atoms.size() * atoms.size(), 99999999.0 );
    for ( size_t row {}; row + 1 < atoms.size(); ++row ) {
      for ( size_t col {row + 1}; col < atoms.size(); ++col ) {
        if ( sametype( atoms[row], atoms[col], similar, chargediff ) ) {
          distmat[ row * atoms.size() + col ] = distance( atoms[row], atoms[col] );
        }
      }
    }
    auto nearest = std::min_element( begin(distmat), end(distmat) );
    best = *nearest;
    auto pos = std::distance( begin(distmat), nearest );
    double limit = cutoff;
    auto it = cutmap.find( atoms[ pos / atoms.size() ].type );
    if ( it != cutmap.end() ) {
      const double limit2 = it.value().toDouble();
      if ( limit2 < limit ) {
        limit = limit2;
      }
      limit = it.value().toDouble();
      }
    // use the smaller cutoff when atomtypes can differ
    it = cutmap.find( atoms[ pos % atoms.size() ].type );
    if ( it != cutmap.end() ) {
      const double limit2 = it.value().toDouble();
      if ( limit2 < limit ) {
        limit = limit2;
      }
    }
    // only atoms within cutoff limit can be merged
    if ( limit < best ) break;

    auto& to = atoms[ pos / atoms.size() ].posit;
    const auto& from = atoms[ pos % atoms.size() ].posit;
    to.insert( to.end(), from.begin(), from.end() );

    // keep the most extreme charge in the cluster
    double& toc = atoms[ pos / atoms.size() ].charge;
    const double frc = atoms[ pos % atoms.size() ].charge;
    if ( std::abs(toc) < std::abs(frc) ) toc = frc;

    atoms.erase( begin(atoms) + (pos % atoms.size()) );
  }

  atoms.erase( std::remove_if( atoms.begin(), atoms.end(),
                               [cmin,nibthreshold](const Atom& x)
                               { return x.posit.size() < cmin &&
                                   std::abs(x.charge) <= nibthreshold;
                               }
                 ),
               atoms.end());
  atoms.erase( std::remove_if( atoms.begin(), atoms.end(),
                               [cminchr,nibthreshold](const Atom& x)
                               { return x.posit.size() < cminchr &&
                                   nibthreshold < std::abs(x.charge);
                               }
                 ),
               atoms.end());
}


void internal_method( std::map<int,std::vector<Atom>>& atomcats, size_t molecule,
                      double cutoff, const QCommandLineParser& parser, bool similar,
                      double chargediff, int argc, char *argv[],
                      unsigned long cmin, unsigned long cminchr, double nibthreshold,
                      QMap<QString, QVariant> cutmap )
{
  std::vector<Atom> atoms;
  size_t original_count {};
  for ( auto& cat : atomcats ) {
    auto& acat = cat.second;
    original_count += acat.size();
    internal_merge( acat, cmin, cminchr, nibthreshold, cutoff, similar, chargediff, cutmap );
    atoms.insert( atoms.end(), begin(acat), end(acat) );
  }

  QString prefix = parser.value( "prefix" );
  std::cout << "# Output from overlap " << qPrintable( QCoreApplication::applicationVersion() ) << '\n';
  std::cout << "# Created: " << qPrintable(QDateTime::currentDateTime().toString()) << '\n';
  std::cout << "# Command:";
  for (int a{}; a < argc; ++a ) std::cout << ' ' << argv[a];
  std::cout << '\n';
  if ( atoms.size() == original_count ) {
    std::cerr << "# Note: No atoms were merged due to overlap\n";
    std::cout << "#\n# Note: No atoms were merged due to overlap\n";
  }
  std::cout << '\n';
  header( std::cout, QString("%1%2").arg(prefix).arg(molecule), atoms.size() );
  unsigned long num {0};
  for ( auto atom : atoms ) {
    ++num;
    print( std::cout, atom, num );
    std::cout << '\n';
  }
  std::cout << '\n';
}


//!
//! Bin atoms according to type
//!
std::map<int,std::vector<Atom>> atoms2bins( const std::vector<QStringList>& atoms,
                                            bool usenib, bool nibneutral, double nibthreshold,
                                            const QStringList& deletelist )
{
  std::map<int,std::vector<Atom>> atomcats;
  for ( size_t e = 0; e < atoms.size(); ++e ) {
    const auto& atom = atoms[e];
    if ( atom.size() == 9 ) {
      QString type = atom[5];
      double charge = atom[8].toDouble();
      if ( usenib ) {
        if ( charge < -nibthreshold ) {
          type = "O.3";
        }
        else if  ( charge > nibthreshold ) {
          type = "N.3";
        }
        else {
          if ( nibneutral ) charge = 0.0;
          if ( type != "C.ar" ){
            type = "C.3";
          }
        }
      }

      if ( ! deletelist.contains( type ) ) {
        int tcat = atomtype(type);
        auto acount = atomcats[tcat].size();
        atomcats[tcat].emplace_back( atom[0],
                                     QString( "%1_%2_%3" ).arg(tcat).arg(acount).arg(atom[1]),
                                     Point{atom[2].toDouble(), atom[3].toDouble(), atom[4].toDouble()},
                                     type, charge );
      }
    }
  }
  return atomcats;
}


//!
//! Output pairs of atoms with similarity (computed from distance with cutoffs)
//!
void bins2abc( std::ostream& ostr, const std::map<int,std::vector<Atom>>& atomcats,
               double cutoff, bool similar,
               double chargediff, QMap<QString, QVariant> cutmap )
{
  for ( const auto& cat : atomcats ) {
    const auto& acat = cat.second;
    for ( size_t row {}; row + 1 < acat.size(); ++row ) {
      double maxdist = cutoff * cutoff;
      auto it = cutmap.find(acat[row].type);
      if ( it != cutmap.end() ) {
        maxdist = it.value().toDouble();
        maxdist *= maxdist;
      }
      for ( size_t col {row + 1}; col < acat.size(); ++col ) {
        if ( sametype( acat[row], acat[col], similar, chargediff ) ) {
          auto d = maxdist - sdist( acat[row], acat[col] );
          if ( 0 < d ) {
            ostr << qPrintable( QString( "%1 %2 %3\n" )
                                .arg(acat[row].name )
                                .arg(acat[col].name )
                                .arg( d ) );
          }
        }
      }
    }
  }
}


int main( int argc, char *argv[] )
{
  QCoreApplication app(argc, argv);
  QCoreApplication::setOrganizationName("SBL");
  QCoreApplication::setApplicationName("o-lap");
  QCoreApplication::setApplicationVersion("2023-08-10");

  QCommandLineParser parser;
  parser.setApplicationDescription("Remove overlapping atoms from a model.\n\n"
                                   "May use Markov Cluster Algorithm (MCL) tool for clustering.\n"
    "Output is either a mol2 model, input for MCL, or atom types in MCL clusters.");
  const auto helpOption = parser.addHelpOption();
  parser.addVersionOption();
  parser.addOption( {"cutoffs", "JSON formatted cutoffs for atom types. Defaults are are from file 'cutoffs.json'.", "file/json" } );
  parser.addOption( {{"c", "cutoff"}, "Cutoff distance. Effective only when shorter than default atomtype specific values (default: 1.1)", "num", "1.1"} );
  parser.addOption( {"showcutoffs", "Show JSON formatted cutoffs for atom types and exit." } );
  parser.addOption( {"similarjson", "JSON formatted atom types.", "file" } );
  parser.addOption( {"showsimilar", "Show similar atom types and exit." } );
  parser.addOption( {{"s", "similar"}, "Cluster similar atom types. Types are similar, if in same category. See '--showsimilar'"} );
  parser.addOption( {"chargediff", "Charges must be within <num> to cluster (default: 0.2).", "num", "0.2"} );
  parser.addOption( {"deletetypes", "Comma-separated list of atom types to discard completely.", "str",} );
  parser.addOption( {"nib", "Convert atom types to (positive) N.3, (negative) O.3, and (neutral) C.3/C.ar based on charge in NIB-like manner."} );
  parser.addOption( {"nibthreshold", "Threshold of charge to bin atoms into N, C, O classes (default: 0.2).", "num", "0.2"} );
  parser.addOption( {"nibcharged", "Keep all charges in model processed with nib option."} );
  parser.addOption( {"clustermin", "Minimum size of cluster to include (default: 1).", "int", "1"} );
  parser.addOption( {"clusterminchr", "Minimum size of cluster for charged atoms. Atom is charged, if abs(charge) exceeds nibthreshold. (default: clustermin)", "int"} );
  parser.addOption( {"abcout", "Create ABC-format input for MCL and exit."} );
  parser.addOption( {"mcl", "Create ABC-format input for MCL and run MCL."} );
  parser.addOption( {"mclI", "MCL main inflation value.", "num"} );
  parser.addOption( {"mclte", "MCL expansion thread number.", "int"} );
  parser.addOption( {"mapmcl", "Map mcl clusters to atoms. The <file> must be output from MCL that corresponds to the model.", "file"} );
  parser.addOption( {"mcltype", "Show types of clustered atoms.  Requires mapmcl."} );
  parser.addOption( {"prefix", "Prefix of the output molecule's name (default: model).", "str", "model"} );
  parser.addPositionalArgument("model", QCoreApplication::translate("main", "Mol2-file"));

  parser.process( app );
  double cutoff = parser.value( "cutoff" ).toDouble();
  double chargediff = parser.value( "chargediff" ).toDouble();
  double nibthreshold = parser.value( "nibthreshold" ).toDouble();
  unsigned long cmin = parser.value( "clustermin" ).toULong();
  unsigned long cminchr = cmin;
  if ( parser.isSet( "clusterminchr" ) ){
    cminchr = parser.value( "clusterminchr" ).toULong();
  }

  bool similar = false;
  if ( parser.isSet( "similar" ) ){
    similar = true;
  }

  QString prefix = parser.value( "prefix" );

  bool usenib = false;
  if ( parser.isSet( "nib" ) ){
    usenib = true;
  }
  bool nibneutral = true;
  if ( parser.isSet( "nibcharged" ) ){
    nibneutral = false;
  }


  QMap<QString, QVariant> cutmap = readjson( "cutoffs.json", parser.value( "cutoffs" ) );
  auto it = cutmap.find("*");
  if ( it != cutmap.end() ) cutoff = it.value().toDouble();

  if ( parser.isSet( "showcutoffs" ) ){
    if ( it == cutmap.end() ) cutmap.insert("*", cutoff);
    for ( const auto& type : atomtypes ) {
      it = cutmap.find(type.first);
      if ( it == cutmap.end() ) cutmap.insert(type.first, cutoff);
    }

    auto obj = QJsonObject::fromVariantMap( cutmap );
    QJsonDocument doc { obj };
    auto arr = doc.toJson( QJsonDocument::Indented );
    std::cout << qPrintable( arr );
    return 0;
  }


  QString similarjson = parser.value( "similarjson" );
  if ( similarjson.isEmpty() )
  {
    similarjson = QStandardPaths::locate( QStandardPaths::AppDataLocation, "atomtypes.json" );
    if ( similarjson.isEmpty() ) {
      QDir myloc( QCoreApplication::applicationDirPath() );
      myloc.cdUp();
      myloc.cd( "share" );
      myloc.cd( QCoreApplication::organizationName() );
      myloc.cd( QCoreApplication::applicationName() );
      if ( myloc.exists() ) {
        similarjson = myloc.path() + "/atomtypes.json";
      }
    }
  }

  if ( ! similarjson.isEmpty() )
  {
    QByteArray ba;
    QFileInfo fi(similarjson);
    if ( fi.isFile() ) {
      QFile cfile(similarjson);
      if (!cfile.open(QIODevice::ReadOnly | QIODevice::Text))
        return 1;
      ba = cfile.readAll();
    }
    QJsonParseError err;
    auto doc = QJsonDocument::fromJson( ba, &err );
    if ( err.error == QJsonParseError::NoError ) {
      atomtypes.clear();
      auto smap = doc.object().toVariantMap();
      for ( auto o = std::begin(smap); o != std::end(smap); ++o ) {
        atomtypes[ qPrintable(o.key()) ] = o.value().toInt();
      }
    }
    else{
      std::cout << "JSON state:" << qPrintable( err.errorString() ) << '\n';
    }
  }

  if ( parser.isSet( "showsimilar" ) ){
    showsimilar();
    return 0;
  }


  QString deletetypes = parser.value( "deletetypes" );
  QStringList deletelist;
  if ( ! deletetypes.isEmpty() ) {
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    deletelist = deletetypes.split(',', QString::SkipEmptyParts);
#else
    deletelist = deletetypes.split(",", Qt::SkipEmptyParts);
#endif
  }


  if ( parser.isSet( "mcltype" ) && ! parser.isSet( "mapmcl" ) ) {
    std::cerr << "Option --mcltype requires --mapmcl.\n";
    return 2;
  }

  const auto positionalArguments = parser.positionalArguments();
  if ( positionalArguments.size() != 1 ) {
    parser.showHelp( 1 );
  }
  else {
    QFile file( positionalArguments.at(0) );
    if ( !file.exists() ) {
      std::cerr << "File " << qPrintable(file.fileName()) << " does not exist.\n";
      return 4;
    }
    if ( !file.open(QIODevice::ReadOnly | QIODevice::Text) ) {
      std::cerr << "Can't open file " << qPrintable(file.fileName()) << "\n";
      return 5;
    }

    QTextStream in( &file );
    auto mols = parse( in );

    for ( size_t i=0; i < mols.size(); ++i )
    {
      auto bins = atoms2bins( mols[i].atoms, usenib, nibneutral,
                              nibthreshold, deletelist );

      QString mcldata = parser.value( "mapmcl" );
      if ( ! mcldata.isEmpty() )
      {
        QFile data( mcldata );
        if ( data.open(QFile::ReadOnly) ) {
          QTextStream input( &data );
          if ( parser.isSet( "mcltype" ) ) {
            QString str;
            QTextStream output( &str );
            mcl2types( input, bins, output );
            std::cout << qPrintable( str );
          }
          else {
            mcl2atoms( input, bins, i, prefix, argc, argv, cmin, cminchr, nibthreshold );
          }
        }
      }
      else if ( parser.isSet( "abcout" ) )
      {
        bins2abc( std::cout, bins, cutoff, similar, chargediff, cutmap );
      }
      else if ( parser.isSet( "mcl" ) )
      {
        std::ostringstream ostr;
        bins2abc( ostr, bins, cutoff, similar, chargediff, cutmap );

        if ( ostr.str().empty() ) {
          std::cerr << "# Note: No atoms were merged due to overlap\n";
        } else {
          // Use 'mcl' for clustering and map result back to atoms
          QProcess mcl;
          QStringList mclopt {"-", "--abc", "-V", "all" };
          if ( parser.isSet( "mclI" ) ) {
            mclopt << "-I" << parser.value( "mclI" );
          }
          if ( parser.isSet( "mclte" ) ) {
            mclopt << "--te" << parser.value( "mclte" );
          }
          mclopt << "-o" << "-";

          mcl.start("mcl", mclopt );
          if ( !mcl.waitForStarted() ) {
            std::cerr << "Failed to start mcl\n";
            return 1;
          }
          mcl.write( ostr.str().c_str() );
          mcl.closeWriteChannel();
          if ( !mcl.waitForFinished( -1 ) ) {
            return 2;
          }
          QTextStream input( mcl.readAllStandardOutput() );
          mcl2atoms( input, bins, i, prefix, argc, argv, cmin, cminchr, nibthreshold );
        }
      }
      else
      {
        // Use iterative "combine nearest" to merge atoms that are within cutoff
        internal_method( bins, i, cutoff, parser, similar, chargediff,
                         argc, argv, cmin, cminchr, nibthreshold, cutmap );
      }
    }
  }
}
