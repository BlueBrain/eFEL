""" eFEL setup """

# pylint: disable=C0325

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

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension  # pylint: disable=E0611,F0401

import os
import versioneer

BASEDIR = os.path.dirname(os.path.abspath(__file__))

cppcore_dir = os.path.join('efel', 'cppcore')
cppcore_sources = ['cppcore.cpp',
                   'Utils.cpp',
                   'LibV1.cpp',
                   'LibV2.cpp',
                   'LibV3.cpp',
                   'LibV4.cpp',
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
                   'LibV4.h',
                   'LibV5.h',
                   'FillFptrTable.h',
                   'DependencyTree.h',
                   'efel.h',
                   'cfeature.h',
                   'Global.h',
                   'mapoperations.h',
                   'types.h',
                   'eFELLogger.h']
cppcore_sources = [
    os.path.join(
        cppcore_dir,
        filename) for filename in cppcore_sources]
cppcore_headers = [
    os.path.join(
        'cppcore',
        filename) for filename in cppcore_headers]

cppcore = Extension('efel.cppcore',
                    sources=cppcore_sources,
                    include_dirs=['efel/cppcore/'])
setup(
    name="efel",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=['numpy>=1.6', 'six'],
    extras_require={'neo': ['neo[neomatlabio]>=0.5.1']},
    packages=['efel'],
    author="BlueBrain Project, EPFL",
    maintainer="Werner Van Geit",
    maintainer_email="werner.vangeit@epfl.ch",
    description="Electrophys Feature Extract Library (eFEL)",
    long_description="The Electrophys Feature Extract Library (eFEL) allows "
    "neuroscientists to automatically extract features from time series data "
    "recorded from neurons (both in vitro and in silico). "
    "Examples are the action potential width and amplitude in "
    "voltage traces recorded during whole-cell patch clamp experiments. "
    "The user of the library provides a set of traces and selects the "
    "features to be calculated. The library will then extract the requested "
    "features and return the values to the user.",
    license="LGPLv3",
    keywords=(
        'feature',
        'extraction',
        'electrophysiology',
        'BlueBrainProject'),
    url='https://github.com/BlueBrain/eFEL',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: GNU Lesser General Public '
        'License v3 (LGPLv3)',
        'Programming Language :: Python :: 2.7',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering',
        'Topic :: Utilities'],
    package_data={
        '': [
            'DependencyV5.txt',
            'VERSION.txt',
            'README.md',
            'GITHASH.txt'] + cppcore_headers},
    ext_modules=[cppcore])
