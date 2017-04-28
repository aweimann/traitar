import os
from setuptools import find_packages
from setuptools import setup
import re


VERSIONFILE=os.path.join('traitar', '_version.py')
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else: 
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

long_description = open('README.rst', 'r').read()

setup(name='traitar',
        version = verstr,
        description='traitar - The microbial trait analyzer',
        long_description = long_description,
        url = 'http://github.com/aweimann/traitar',
        author='Aaron Weimann',
        author_email='weimann@hhu.de',
        license='GNU General Public License, version 3 (GPL-3.0)',
        packages= ['traitar'],
        include_package_data = True,
        scripts = ['bin/traitar', 'bin/merge_preds', 'bin/heatmap', 'bin/domtblout2gene_generic', 'bin/predict', 'bin/hmmer2filtered_best', 'bin/hmm2gff'],
        zip_safe=False,
        install_requires = ["pandas >= 0.13.1", "matplotlib >= 1.3.1", "numpy >= 1.6", "scipy >= 0.13.3"])
