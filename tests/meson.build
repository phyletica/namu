catch2_proj = subproject('catch2')

catch2_with_main_dep = catch2_proj.get_variable('catch2_with_main_dep')

phycoeval_test_files = files(
  'phycoeval_testing.cpp',
  'test_tree.cpp',
  'test_node.cpp',
)

phycoeval_test_exe = executable(
    'test_phycoeval',
    phycoeval_test_files,
    dependencies : [ phycoeval_dep, catch2_with_main_dep ]
)

test('phycoeval test', phycoeval_test_exe)
