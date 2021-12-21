import os
import platform
import setuptools
import sys
import subprocess
import shutil
from pathlib import Path
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

from distutils.util import get_platform

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir = ""):
        Extension.__init__(self, name, sources = [])
        self.sourcedir = Path(sourcedir).absolute()


class CMakeBuild(build_ext):
    """
    This class is built upon https://github.com/diegoferigo/cmake-build-extension/blob/master/src/cmake_build_extension/build_extension.py and https://github.com/pybind/cmake_example/blob/master/setup.py
    """

    user_options = build_ext.user_options + [
        ("define=", "D", "Define variables for CMake")
    ]

    def initialize_options(self):
        super().initialize_options()
        self.define = None

    def finalize_options(self):
        # Parse the custom CMake options and store them in a new attribute
        defines = [] if self.define is None else self.define.split(";")
        self.cmake_defines = [f"-D{define}" for define in defines]

        super().finalize_options()

    def build_extension(self, ext: CMakeExtension):
        extdir = str(Path(self.get_ext_fullpath(ext.name)).parent.absolute())

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
        ninja_path = str(shutil.which('ninja'))

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        configure_args = [
            "-GNinja",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            f"-DCMAKE_MAKE_PROGRAM={ninja_path}",
        ]
        configure_args += self.cmake_defines
        
        build_args = []


        # Add more platform dependent options
        if platform.system() == "Darwin":
            pass
        elif platform.system() == "Linux":
            if platform.machine() == "x86_64":
                configure_args += [
                    "-DENABLE_AVX=ON"
                ] # Enable AVX if x64 on Linux
        elif platform.system() == "Windows":
            pass
        else:
            raise RuntimeError(f"Unsupported '{platform.system()}' platform")

        if not Path(self.build_temp).exists():
            os.makedirs(self.build_temp)
        
        subprocess.check_call(
            ["cmake", str(ext.sourcedir)] + configure_args, cwd = self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd = self.build_temp
        )


with open("union_find/_version.py") as f:
    version = f.readlines()[-1].split()[-1].strip("\"'")

requirements = [
    "ninja"
]

info = {
    "name": "UnionFindCPP",
    "version": version,
    "maintainer": "Chae-Yeun Park",
    "maintainer_email": "chae.yeun.park@gmail.com",
    "url": "https://github.com/chaeyeunpark/UnionFind",
    "license": "GNU Lesser General Public License v2.1",
    "packages": find_packages(where="."),
    "package_data": {"union_find": ["cpp/*"]},
    "entry_points": {
    },
    "description": "C++ implemenation of UnionFind with a python binding",
    "long_description": open("README.md").read(),
    "long_description_content_type": "text/markdown",
    "provides": ["union_find"],
    "install_requires": requirements,
    "ext_modules": [CMakeExtension("_union_find_py")],
    "cmdclass": {"build_ext": CMakeBuild},
    "ext_package": "union_find",
}

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Operating System :: POSIX :: Linux",
    "Operating System :: Microsoft :: Windows",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering :: Physics",
]

setup(classifiers=classifiers, **(info))
