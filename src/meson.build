executable('namu', 'main.cpp',
    link_with : namu_lib,
    dependencies : [ ncl_lib, boost_program_options_lib, boost_filesystem_lib, beagle_lib ],
    include_directories : [ namu_include, boost_include, ncl_include, beagle_include ],
    install : true,
    install_dir : './bin')
