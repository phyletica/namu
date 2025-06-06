project('phycoeval', 'cpp',
  default_options : ['cpp_std=c++14', 'cpp_link_args=-static', 'prefix=' + meson.current_source_dir() + '/distr'],
  version : '0.1')

cpp = meson.get_compiler('cpp')

ncl_lib = cpp.find_library('ncl',
    dirs: [meson.current_source_dir() + '/dependencies/ncl-build/installed/lib/static'],
    required: true
)
ncl_include = include_directories('./dependencies/ncl-build/installed/include')

boost_include = include_directories('./dependencies/boost_1_71_0', is_system : true)

phycoeval_include = include_directories('include')

phycoeval_lib = static_library(
    'phycoeval_lib',
    'src/main.cpp',
    include_directories : [ phycoeval_include, boost_include, ncl_include ],
    install : true
)

phycoeval_dep = declare_dependency(
    link_with : phycoeval_lib,
    include_directories : phycoeval_include
)

subdir('tests')
subdir('src')
