""" eFEL setup """

"""
Copyright (c) 2015, EPFL/Blue Brain Project

 This file is part of eFEL <https://github.com/BlueBrain/eFEL>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


from setuptools import setup, Extension

import os
import versioneer

BASEDIR = os.path.dirname(os.path.abspath(__file__))

cppcore_dir = os.path.join('efel', 'cppcore')
cppcore_sources = ['cppcore.cpp',
                   'Utils.cpp',
                   'LibV1.cpp',
                   'LibV2.cpp',
                   'LibV3.cpp',
                   'LibV5.cpp',
                   'FillFptrTable.cpp',
                   'DependencyTree.cpp',
                   'efel.cpp',
                   'cfeature.cpp',
                   'mapoperations.cpp']
cppcore_headers = ['Utils.h',
                   'LibV1.h',
                   'LibV2.h',
                   'LibV3.h',
                   'LibV5.h',
                   'FillFptrTable.h',
                   'DependencyTree.h',
                   'efel.h',
                   'cfeature.h',
                   'Global.h',
                   'mapoperations.h',
                   'types.h',
                   'eFELLogger.h',
                   'EfelExceptions.h']
cppcore_sources = [
    os.path.join(
        cppcore_dir,
        filename) for filename in cppcore_sources]
cppcore_headers = [
    os.path.join(
        'cppcore',
        filename) for filename in cppcore_headers]


coverage_flags = []
if os.environ.get('EFEL_COVERAGE_BUILD'):
    coverage_flags = ['-fprofile-arcs', '-ftest-coverage']

cppcore = Extension('efel.cppcore',
                    sources=cppcore_sources,
                    include_dirs=['efel/cppcore/'],
                    extra_compile_args=coverage_flags + ['-std=c++17'],
                    extra_link_args=coverage_flags)

with open("README.md", encoding="utf-8") as f:
    README = f.read()

setup(
    name="efel",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=['numpy>=1.6', 'neo>=0.5.2'],
    packages=['efel', 'efel.pyfeatures', 'efel.units'],
    author="BlueBrain Project, EPFL",
    maintainer="Werner Van Geit",
    maintainer_email="werner.vangeit@epfl.ch",
    description="Electrophys Feature Extract Library (eFEL)",
    long_description=README,
    long_description_content_type="text/markdown",
    license="LGPLv3",
    keywords=[
        'feature',
        'extraction',
        'electrophysiology',
        'BlueBrainProject'],
    url='https://github.com/BlueBrain/eFEL',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: GNU Lesser General Public '
        'License v3 (LGPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering',
        'Topic :: Utilities'],
    package_data={
        '': ['DependencyV5.txt',
             'README.md',
             'units/units.json'] + cppcore_headers},
    ext_modules=[cppcore])
