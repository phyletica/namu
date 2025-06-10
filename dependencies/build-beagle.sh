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

beagle_name="beagle-lib-3.1.2"
beagle_dir="${dep_dir}/${beagle_name}"
install_dir="${beagle_dir}/installed"

if [ -d "$beagle_dir" ]
then 
    echo "ERROR: beagle directory '$beagle_dir' already exists."
    echo "To replace, please remove this directory and re-run this script."
    exit 1
fi

if [ ! -d "$dep_dir" ]
then
    mkdir "$dep_dir"
fi

cd "$dep_dir"

echo
echo "Downloading beagle release 3.1.2..."
echo
curl -L --output "${beagle_name}.tar.gz" "https://github.com/beagle-dev/beagle-lib/archive/v3.1.2.tar.gz"
echo
echo "Extracting beagle..."
echo
tar xzf "${beagle_name}.tar.gz"

echo
echo "Compiling beagle..."
echo
mkdir "$install_dir"
cd "$beagle_dir"
./autogen.sh
./configure --prefix="$install_dir" --with-jdk=no CXXFLAGS="-std=c++11"
make
make install

env_path="${dep_dir}/beagle-env.sh"
echo export BEAGLE_PREFIX="${install_dir}" > "$env_path"
echo export LD_LIBRARY_PATH="${install_dir}/lib:\${LD_LIBRARY_PATH}" >> "$env_path"
