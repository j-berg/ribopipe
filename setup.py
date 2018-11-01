"""
RiboPipe
An assembly and analysis pipeline for sequencing data
alias: ripopipe 

Copyright (C) 2018  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from setuptools import setup
import re

with open('ribopipe/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)



setup(
    name = 'RiboPipe',
    version = version,
    description = 'A Ribosome Profiling and Small RNA Data Handling Pipeline',
    long_description = open('README.md').read(),
    author = 'Jordan Berg',
    author_email = 'jordan.berg@biochem.utah.edu',
    url = 'https://github.com/j-berg/ribopipe',
    packages = ['ribopipe'],
    package_dir = {'ribopipe': 'ribopipe'},
    include_package_data = True,
    license = 'GPL-3.0',
    zip_safe = False,

    entry_points = {
        'console_scripts': [
            'ribopipe = ribopipe.__main__:main'
        ]
    },

    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)
