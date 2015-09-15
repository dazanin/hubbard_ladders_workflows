#!/bin/bash

# Dataverse data
# http://dx.doi.org/10.7910/DVN/I5ANSU
URL=https://dataverse.harvard.edu/api/access/datafile/2704552

# Archive disk
# URL=http://archive.comp-phys.org/phys.ethz.ch/dolfim/Papers/hubbard_ladder/data_extracted.tar.gz

# Checksum
FILECHECK=fd178d8f27c8b8339d22016fdfcec8f2
# Filename
ONAME=data_extracted.tar.gz
DIST=data_extracted

if [ "$(ls -A $DIST)" ]; then
     echo "The $DIST directory is NOT empty."
     read -p "Do you want to continue? [Yn] " prompt
     if [[ $prompt =~ [nN](o)* ]]; then
         exit
     fi
fi

echo "Downloading file..."
curl -o $ONAME $URL

echo "Verifying file..."
MYCHECK=`openssl md5 -r $ONAME | awk '{print $1}'`
if [ $MYCHECK != $FILECHECK ]; then
    echo "Error: the downloaded file does not match the expected checksum."
    echo "To continue modify the script or contact the package maintainer."
    exit
fi


echo "Uncompressing data..."
tar xzf $ONAME

echo "Done."

