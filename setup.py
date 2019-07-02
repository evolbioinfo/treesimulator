import os
from setuptools import setup, find_packages

setup(
    name='treesimulator',
    packages=find_packages(),
    include_package_data=True,
    package_data={'treesimulator': [os.path.join('..', 'README.md')]},
    long_description=open('README.md').read(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.0.1',
    description='Simulation of rooted phylogenetic trees following an epidemiological model.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://gitlab.pasteur.fr/phylo/treesimulator',
    keywords=['phylogeny', 'simulation'],
    install_requires=['ete3', 'pandas', 'numpy', 'scipy'],
    entry_points={
            'console_scripts': [
                'treesimulator = treesimulator.main:main',
            ]
    },
)
