#!/usr/bin/env bash

archive_dir="./archives"
if [[ ! -d ./archives ]] ; then
  mkdir $archive_dir
fi
# exit if sphinx-build not installed
if ! command -v sphinx-build &> /dev/null
then
    echo "sphinx-build could not be found"
    echo "might you havn't entered the virtual environment"
    exit 1
fi
# build and upload the archive
cd docs/
make html
cd build

# Create tar.gz archive of the html folder with a timestamp
timestamp=$(date +%Y%m%d%H%M%S)
archive_name="$archive_dir/html_docs_$timestamp.tar.gz"
tar -czf $archive_name html

# Create docs folder on the remote host if it doesn't exist
# requires the /opt/statics to be writeable by the ssh user
if [[ -z $1 ]] ; then
  echo "No remote host provided, only archive created"
  echo "upload skipped"
  exit 0
fi

echo "Archive: $archive_name"
ssh $1 "mkdir -p /opt/statics/docs"

# Upload the archive to the remote host using rsync
rsync -avz $archive_name $1:/opt/statics/docs

# Extract the archive on the remote host
ssh $1 "tar -xzf /opt/statics/docs/$archive_name -C /opt/statics/docs"
