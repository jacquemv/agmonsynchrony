import platform
from setuptools import setup, Extension
from Cython.Distutils import build_ext

NAME = "agmonsynchrony"
VERSION = "0.1.0"
DESCR = "Synchrony index between time series based on Agmon's paper"
KEYWORDS = "time series,synchrony,neuron,spikes"
URL = "http://github.com/jacquemv/agmonsynchrony"
REQUIRES = ['numpy', 'cython']
AUTHOR = "Vincent Jacquemet"
EMAIL = "vincent.jacquemet@umontreal.ca"
LICENSE = "MIT"
SRC_DIR = "agmonsynchrony"
PACKAGES = [SRC_DIR]

if platform.system() == 'Windows':
    compiler_args = ['/O2']
    linker_args = []
else:
    compiler_args = ['-O3']
    linker_args = []
    
ext = Extension(SRC_DIR + ".synchrony",
                sources=[SRC_DIR + "/windowoverlap.cpp",
                         SRC_DIR + "/synchrony.pyx"],
                libraries=[],
                extra_compile_args=compiler_args,
                extra_link_args=linker_args,
                language="c++",
                include_dirs=[SRC_DIR])
ext.cython_directives = {'language_level': "3"}

EXTENSIONS = [ext]

setup(install_requires=REQUIRES,
      packages=PACKAGES,
      zip_safe=False,
      name=NAME,
      version=VERSION,
      description=DESCR,
      keywords=KEYWORDS,
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      author=AUTHOR,
      author_email=EMAIL,
      url=URL,
      license=LICENSE,
      cmdclass={"build_ext": build_ext},
      ext_modules=EXTENSIONS
)
