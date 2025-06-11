catch2_proj = subproject('catch2', default_options : ['tests=false'])

catch2_with_main_dep = catch2_proj.get_variable('catch2_with_main_dep')

namu_test_files = files(
  'namu_testing.cpp',
  'test_tree.cpp',
  'test_node.cpp',
  'test_string_util.cpp',
  'test_error.cpp',
)

namu_test_exe = executable(
    'test_namu',
    namu_test_files,
    link_with : namu_lib,
    dependencies : [ catch2_with_main_dep, ncl_lib, boost_program_options_lib, boost_filesystem_lib, beagle_lib ],
    include_directories : [ namu_include, boost_include, ncl_include, beagle_include, eigen_include ]
)

test('namu test', namu_test_exe)
