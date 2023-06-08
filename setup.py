import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 0
REV = 1


VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'GWAS/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='GWAS-py',
    version=VERSION,
    description='GWAS-py Project: a linear GWAS tool with plotting utility',
    author='Cathy(Weiwen), Jenelle, Jiayi',
    author_email='w3dong@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "GWAS-py=GWAS.GWAS:main"
        ],
    },
)