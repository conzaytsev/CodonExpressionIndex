from setuptools import setup


setup(
    name="cei", 
    version='0.1.1',
    author="Konstantin Zaytsev",
    url='https://github.com/conzaytsev/CodonExpressionIndex',
    description='Codon Expression Index',
    long_description='Python module for analysis of codon influence on protein expression.',
    packages=['cei'],
    install_requires=['scipy'],
    python_requires='>3.9',
    package_data={'cei': ['datasets/*.csv']},
    classifiers= [
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
    ]
)
