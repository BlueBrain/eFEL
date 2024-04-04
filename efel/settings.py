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

from pathlib import Path


class Settings:
    """FEL settings class"""

    def __init__(self):
        self.threshold = -20.0
        self.derivative_threshold = 10.0
        self.down_derivative_threshold = -12.0
        self.dependencyfile_path = str(
            Path(__file__).parent.absolute() / "DependencyV5.txt"
        )
