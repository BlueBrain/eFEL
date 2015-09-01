"""eFEL version"""

# pylint: disable=W0633, W0702

import os
import inspect

version_filename = 'VERSION.txt'
githash_filename = 'GITHASH.txt'


def _get_git_version(script_path):
    """Get version from git"""

    # Make sure we're not a git repository on much higher level
    if os.path.exists(os.path.join(script_path, '../.git')):
        try:
            import subprocess
            p = subprocess.Popen(
                ['git', 'describe', '--tags', '--long'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=script_path)
            return p.stdout.readlines()[0].split('-')
        except:
            return None
    else:
        return None


def _get_version_number():
    """Get the current version"""

    script_filename = inspect.getframeinfo(inspect.currentframe()).filename
    script_path = os.path.dirname(os.path.abspath(script_filename))
    version_fullpath = os.path.join(script_path, version_filename)
    githash_fullpath = os.path.join(script_path, githash_filename)

    git_version_return = _get_git_version(script_path)
    if git_version_return is None:
        if os.path.exists(version_fullpath):
            with open(version_fullpath, 'r') as version_file:
                version_number = version_file.readline()
        else:
            print("WARNING: No valid version number found for efel")
            return None
    else:
        major_version, rev_index, githash = git_version_return
        version_number = '%s.%s' % (major_version, rev_index)
        with open(version_fullpath, 'w') as version_file:
            version_file.write(version_number)
        with open(githash_fullpath, 'w') as githash_file:
            githash_file.write(githash)

    return version_number

version = _get_version_number()
__version__ = version
VERSION = version
