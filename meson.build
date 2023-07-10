project('mosim', 'c', 'cpp', 
  default_options : ['cpp_std=c++20'], 
  version: '1.0'
)
cmake = import('cmake')

if get_option('buildtype').startswith('debug')
  add_project_arguments('-DDEBUG', language : 'cpp')
endif

fmt_dep = dependency('fmt')
eigen_dep = dependency('eigen3')
glfw_dep = dependency('glfw3')
glm_dep = dependency('glm')
glew_dep = dependency('glew')
omp_dep = dependency('openmp')
libarchive_dep = dependency('libarchive')

sources = [
  'src/main.cpp',
  'src/hf/hf.cpp',
  'src/hf/basis.cpp',
  'src/init.cpp',
  'src/gl/orbital_render.cpp',
  'src/gl/render.cpp',
  'src/gl/shader.cpp',
  'src/gl/window.cpp',
  'src/resources.cpp',
  'src/hf/chem.cpp',
  'src/gl/font.cpp'
]

include_dirs = [
  'include',
  '.'
]

tar =  find_program('tar')
xxd = find_program('xxd')

asset_dir = meson.source_root()/'assets'

assets_compressed = custom_target('assets-compressed',
  output: 'assets.tar',
  input: [ ],
  command : [tar, '-cvf', '@OUTPUT@', '-C', meson.source_root(), 'assets'],
  build_always: true,
  build_by_default: true
)

assets_header = custom_target('assets-header',
  output: 'assets.h',
  input: assets_compressed,
  command : [xxd, '-i', '@INPUT@', '@OUTPUT@'],
  build_always: true,
  build_by_default: true
)

stacktrace_opts = cmake.subproject_options()
stacktrace_opts.add_cmake_defines({'USE_STACKTRACE_BACKEND': 'execinfo', 'USE_DECODE_BACKEND': 'libbfd'})
stacktrace_subproj = cmake.subproject('stacktrace', options: stacktrace_opts)
stacktrace_dep = stacktrace_subproj.dependency('stacktrace')

nfd_subproj = cmake.subproject('nativefiledialog-extended')
nfd_dep = nfd_subproj.dependency('nfd')

conf_data = configuration_data()
conf_data.set_quoted('VERSION', meson.project_version())
conf_data.set_quoted('PROJECT_NAME', 'mosim')
conf_data.set_quoted('MESON_C_COMPILER', meson.get_compiler('c').get_id())
conf_data.set_quoted('MESON_CXX_COMPILER', meson.get_compiler('cpp').get_id())

configure_file(input: 'build_config.h.in', output: 'build_config.h', configuration : conf_data)

mosim = executable('mosim', sources, 
  dependencies: [fmt_dep, eigen_dep, glfw_dep, glew_dep, glm_dep, stacktrace_dep, libarchive_dep, omp_dep, nfd_dep],
  cpp_pch: 'pch/pch.h', 
  cpp_args : '-DFMT_HEADER_ONLY', 
  link_args: ['-lbfd'],
  include_directories: include_directories(include_dirs),
)

