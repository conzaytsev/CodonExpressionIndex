from setuptools import setup

VERSION = '0.1.0' 
DESCRIPTION = 'Codon Expression Index'
LONG_DESCRIPTION = 'Python module for analysis of codon influence on protein expression.'

# Setting up
setup(
    name="cei", 
    version=VERSION,
    author="Konstantin Zaytsev",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages='cei,
    install_requires=['scipy'],
    classifiers= ["Programming Language :: Python :: 3"]
)
