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

boost_dir="${dep_dir}/${boost_name}"

env_path="${dep_dir}/boost-env.sh"
echo export BOOST_PREFIX="${boost_dir}" > "$env_path"
echo export LD_LIBRARY_PATH="${boost_dir}:\${LD_LIBRARY_PATH}" >> "$env_path"
