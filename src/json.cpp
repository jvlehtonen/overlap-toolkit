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

#include "json.h"


void json2map( QByteArray& ba, QMap<QString, QVariant>& cutmap )
{
  QJsonParseError err;
  auto doc = QJsonDocument::fromJson( ba, &err );
  if ( err.error == QJsonParseError::NoError ) {
    auto obj = doc.object();
    QMap<QString, QVariant> data = obj.toVariantMap();
    auto i = data.constBegin();
    while (i != data.constEnd()) {
      cutmap.insert( i.key(), i.value() );
      ++i;
    }
  }
  else {
    std::cout << "JSON state:" << qPrintable( err.errorString() ) << '\n';
  }
}


//!
//! Read JSON from 'filename' and add userdata.
//! Return two-column table
//!
QMap<QString, QVariant> readjson( const QString& filename, const QString& userdata )
{
  QMap<QString, QVariant> cutmap;
  QByteArray ba;

  // read default cutoffs
  QString appdir = QStandardPaths::locate( QStandardPaths::AppDataLocation, filename );
  if ( appdir.isEmpty() ) {
    //std::cout << "Not in QStandardPaths::AppDataLocation\n";
    QDir myloc( QCoreApplication::applicationDirPath() );
    myloc.cdUp();
    myloc.cd( "share" );
    myloc.cd( QCoreApplication::organizationName() );
    myloc.cd( QCoreApplication::applicationName() );
    if ( myloc.exists() ) {
      //std::cout << "DataLocation: " << qPrintable(myloc.path()) << '\n';
      appdir = myloc.path() + "/" + filename;
    }
  }

  if ( not appdir.isEmpty() )
  {
    QFile cfile(appdir);
    if ( cfile.open(QIODevice::ReadOnly | QIODevice::Text) )
    {
      ba = cfile.readAll();
      json2map( ba, cutmap );
    }
  }

  // add custom cutoffs from user
  if ( ! userdata.isEmpty() )
  {
    QFileInfo fi(userdata);
    if ( fi.isFile() ) {
      QFile cfile(userdata);
      if (!cfile.open(QIODevice::ReadOnly | QIODevice::Text))
        return cutmap;
      ba = cfile.readAll();
    }
    else {
      ba = QByteArray(qPrintable(userdata));
    }
    json2map( ba, cutmap );
  }

  return cutmap;
}
