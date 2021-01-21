"""efel Settings class"""

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

import inspect
import os


def _get_script_path():
    """Get directory path of current script"""
    script_filename = inspect.getframeinfo(inspect.currentframe()).filename
    script_path = os.path.dirname(os.path.abspath(script_filename))

    return script_path


class Settings(object):

    """FEL settings class"""

    def __init__(self):
        self.threshold = -20.0
        self.derivative_threshold = 10.0
        self.down_derivative_threshold = -12.0
        self.dependencyfile_path = os.path.join(
            _get_script_path(),
            'DependencyV5.txt')
