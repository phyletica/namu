project('namu', 'cpp',
  default_options : ['cpp_std=c++14', 'prefix=' + meson.current_source_dir() + '/distr'],
  version : '0.1')

cpp = meson.get_compiler('cpp')

boost_program_options_lib = cpp.find_library('boost_program_options',
    dirs: [meson.current_source_dir() + '/dependencies/boost_1_71_0/stage/lib/static'],
    required: true
)
boost_filesystem_lib = cpp.find_library('boost_filesystem',
    dirs: [meson.current_source_dir() + '/dependencies/boost_1_71_0/stage/lib/static'],
    required: true
)

ncl_lib = cpp.find_library('ncl',
    dirs: [meson.current_source_dir() + '/dependencies/ncl-build/installed/lib/static'],
    required: true
)

beagle_lib = cpp.find_library('hmsbeagle',
    dirs: [meson.current_source_dir() + '/dependencies/beagle-lib-3.1.2/installed/lib'],
    required: true
)

ncl_include = include_directories('./dependencies/ncl-build/installed/include')

boost_include = include_directories('./dependencies/boost_1_71_0', is_system : true)

beagle_include = include_directories('./dependencies/beagle-lib-3.1.2/installed/include/libhmsbeagle-1')

namu_include = include_directories('include')

namu_lib = static_library(
    'namu',
    'src/main.cpp',
    include_directories : [ namu_include, boost_include, ncl_include, beagle_include ],
    install : true
)

# namu_dep = declare_dependency(
#     link_with : namu_lib,
#     include_directories : namu_include
# )

subdir('tests')
subdir('src')
