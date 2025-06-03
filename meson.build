project('phycoeval', 'cpp',
  default_options : ['cpp_std=c++14', 'cpp_link_args=-static', 'prefix=/home/jamie/Dropbox/software/dev/phycoeval/distr'],
  version : '0.1')

phycoeval_include = include_directories('include')

phycoeval_lib = static_library(
    'phycoeval_lib',
    'src/main.cpp',
    include_directories : [ phycoeval_include ],
    install : true
)

phycoeval_dep = declare_dependency(
    link_with : phycoeval_lib,
    include_directories : phycoeval_include
)

subdir('tests')
subdir('src')
