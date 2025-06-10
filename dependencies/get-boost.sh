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

boost_name="boost_1_71_0"
boost_dir="${dep_dir}/${boost_name}"

if [ -d "$boost_dir" ]
then 
    echo "ERROR: boost directory '$boost_dir' already exists."
    echo "To replace, please remove this directory and re-run this script."
    exit 1
fi

if [ ! -d "$dep_dir" ]
then
    mkdir "$dep_dir"
fi

cd "$dep_dir"

echo
echo "Downloading Boost release 1.71.0..."
echo
curl -LO "https://archives.boost.io/release/1.71.0/source/${boost_name}.tar.gz"
echo
echo "Extracting Boost release 1.71.0..."
echo
tar xzf "${boost_name}.tar.gz"

echo
echo "Compiling Boost program_options library..."
echo
cd "$boost_dir"
./bootstrap.sh --with-toolset=gcc --with-libraries=program_options,filesystem,system
./b2 cxxflags="-std=c++11" -d 2

boost_lib_dir="${boost_dir}/stage/lib"
boost_static_dir="${boost_dir}/stage/lib/static"
mkdir "$boost_static_dir"
cp "${boost_lib_dir}"/*.a "$boost_static_dir"

echo
echo "Static libraries copied to:"
echo "$boost_static_dir"
echo

env_path="${dep_dir}/boost-env.sh"
echo export BOOST_PREFIX="${boost_dir}" > "$env_path"
echo export LD_LIBRARY_PATH="${boost_dir}:\${LD_LIBRARY_PATH}" >> "$env_path"
echo export boost_static_dir="${boost_static_dir}" >> "$env_path"
