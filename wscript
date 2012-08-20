VERSION = '0.0.1'
APPNAME = 'wavelet_matrix'

srcdir = '.'
blddir = 'build'

def set_options(ctx):
  ctx.tool_options('compiler_cxx')
  ctx.tool_options('unittestt')

def configure(ctx):
  ctx.check_tool('compiler_cxx')
  ctx.check_tool('unittestt')	
  ctx.check_cxx(lib = 'gtest_main', uselib_store = 'gtest')
  ctx.env.CXXFLAGS += ['-O2', '-Wall', '-W', '-g']

import Scripting
Scripting.dist_exts += ['.sh']

def dist_hook():
    import os
    os.remove('googlecode_upload.py')

def build(bld):
  bld.recurse('src test performance_test')
  bld.install_files('/usr/lib/pkgconfig', 'wavelet_matrix.pc')
