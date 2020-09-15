#This file is forked from https://github.com/pybind/pbtest, original author: Sylvain Corlay

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools
import os, sys

# get __version__, __author__, and __email__
exec(open("./sixdegrees/metadata.py").read())

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

ext_modules = [
    Extension(
        '_sixdegrees',
        [ 
            '_sixdegrees/_sixdegrees.cpp', 
            '_sixdegrees/Utilities.cpp', 
            '_sixdegrees/ssmh.cpp', 
            '_sixdegrees/kleinberg.cpp', 
            '_sixdegrees/original_small_world.cpp', 
            '_sixdegrees/modified_small_world.cpp', 
            '_sixdegrees/random_geometric_small_world.cpp', 
            '_sixdegrees/twoD_random_geometric_kleinberg.cpp', 
            '_sixdegrees/find_first_random_edge.cpp', 
        ],
        include_dirs=[
            get_pybind_include(),
            get_pybind_include(user=True),
            "./_sixdegrees/"
        ],
        language='c++',
    ),
]

def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    fd, fname = tempfile.mkstemp('.cpp', 'main', text=True)
    with os.fdopen(fd, 'w') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
        return False
    return True

def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is preferred over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7','-ftemplate-depth=1024']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(
    name='sixdegrees',
    version=__version__,
    author=__author__,
    author_email=__email__,
    url='https://github.com/benmaier/sixdegrees',
    license=__license__,
    description='Creates modular hierarichical random networks, Kleinberg networks and small world networks in a fast manner.',
    long_description='',
    packages = setuptools.find_packages(),
    ext_modules=ext_modules,
    setup_requires=['pybind11>=2.2.4'],
    install_requires=[
                 'pybind11>=2.2.4',
                 'numpy>=1.15.4',
                 'scipy>=1.2.0',
                 ],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)
