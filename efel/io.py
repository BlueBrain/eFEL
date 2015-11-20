"""IO handler for eFEL"""

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

import os

# Python 2 has urlparse module, Python 3 has urllib.parse
try:
    import urlparse as up
except ImportError:
    # pylint:disable=E0611, F0401
    import urllib.parse as up
    # pylint:enable=E0611,F0401

import mimetypes


def load_fragment(fragment_url, mime_type=None):
    """Load fragment

    Load a fragment (e.g. time series data) from a given URL
    """

    parsed_url = up.urlparse(fragment_url)

    scheme = parsed_url.scheme
    server_loc = parsed_url.netloc
    path = parsed_url.path
    fragment_string = parsed_url.fragment

    if mime_type is None:
        mimetypes.init()
        mime_type, _ = mimetypes.guess_type(path)
        if mime_type is None:
            raise TypeError(
                'load_fragment: impossible to guess MIME type from url, '
                'please specify the type manually as argument: %s' % path)

    if scheme == 'file':
        file_handle = open(os.path.join(server_loc, path), 'r')

    if 'text/' in mime_type:
        import numpy
        if fragment_string == '':
            cols = None
        else:
            import re
            match = re.match("col=([0-9]+)", fragment_string)
            if match is None or len(match.groups()) != 1:
                raise TypeError(
                    "load_fragment: don't understand url fragment %s" %
                    fragment_string)
            else:
                cols = int(match.groups()[0]) - 1

        # Unfortunately we need this if statement
        # Setting usecols to None throws an error in the loadtxt call
        if cols is not None:
            fragment_content = numpy.loadtxt(file_handle, usecols=[cols])
        else:
            fragment_content = numpy.loadtxt(file_handle)

        return fragment_content
    else:
        raise TypeError('load_fragment: unknown mime type %s' % mime_type)
