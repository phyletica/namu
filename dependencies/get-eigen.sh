#!/usr/bin/env bash
 
set -e

# make sure we get back to directory of caller
current_dir="$(pwd)"
function return_on_exit () {
    cd "$current_dir"
}
trap return_on_exit EXIT

# get location of script
dep_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

eigen_version="3.4.0"
eigen_name="eigen-${eigen_version}"
eigen_dir="${dep_dir}/${eigen_name}"

if [ -d "$eigen_dir" ]
then 
    echo "ERROR: eigen directory '$eigen_dir' already exists."
    echo "To replace, please remove this directory and re-run this script."
    exit 1
fi

if [ ! -d "$dep_dir" ]
then
    mkdir "$dep_dir"
fi

cd "$dep_dir"

echo
echo "Downloading eigen release ${eigen_version}..."
echo
curl -LO "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
echo
echo "Extracting eigen release ${eigen_version}..."
echo
tar xzf "${eigen_name}.tar.gz"
