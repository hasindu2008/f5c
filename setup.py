#adapted from https://github.com/lh3/minimap2/blob/master/setup.py

try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension


cmdclass={}

from Cython.Build import build_ext
module_src = 'python/pyf5c.pyx'
cmdclass['build_ext'] = build_ext


sources=[module_src, 'src/align.c', 'src/eventalign.c', 'src/events.c', 'src/f5c.c', 'src/f5cio.c', 'src/freq.c', 'src/freq_merge.c', 'src/hmm.c', 'src/index.c', 'src/meth.c', 'src/meth_main.c', 'src/model.c', 'src/nanopolish_fast5_io.c', 'src/nanopolish_read_db.c', 'src/profiles.c']
depends=['python/pyf5c.pxd', 'python/pyf5c.h', 'src/config.h', 'src/error.h', 'src/f5c.h', 'src/f5cmisc.h', 'src/fast5lite.h', 'src/khash.h', 'src/ksort.h', 'src/logsum.h', 'src/matrix.h', 'src/model.h', 'src/nanopolish_read_db.h', 'src/profiles.h', 'src/str.h']
extra_compile_args = ['-g', '-Wall', '-O2', '-std=c++11', '-Wno-strict-prototypes']
libraries = ['z','m','pthread', 'hdf5_serial', 'hts']
include_dirs = ['htslib']
library_dirs = ['htslib']

#py_inc = [get_python_inc()]

#np_lib = os.path.dirname(numpy.__file__)
#np_inc = [os.path.join(np_lib, 'core/include')]
#cmdclass = {'build_py': build_py}


#cmdclass.update({'build_ext': build_ext})
#packages=['test']


extensions = [Extension('pyf5c',
                  sources = sources,
                  depends = depends,
                  extra_compile_args = extra_compile_args,
                  libraries = libraries,
                  include_dirs = include_dirs,
                  library_dirs = library_dirs,
                  language = 'c++' )]

setup(name = 'f5c',
      version='0.6',
      url = 'https://github.com/hasindu2008/f5c',
      #requires=['numpy (>=1.3.0)'],
      description='f5c python binding',
      author='Hasindu Gamaarachchi',
      author_email='hasindu@garvan.org.au',
      maintainer='Hasindu Gamaarachchi',
      maintainer_email='hasindu@garvan.org.au',
      license = 'MIT',
      keywords = ['nanopore','methylation','event-alignment','signal-alignment','sequence-alignment'],
      #packages=packages,
      cmdclass=cmdclass,
      ext_modules=extensions
      #ext_modules=cythonize(extensions),
      )
