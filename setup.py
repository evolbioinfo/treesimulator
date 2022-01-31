import os
from setuptools import setup, find_packages

setup(
    name='treesimulator',
    packages=find_packages(),
    include_package_data=True,
    package_data={'treesimulator': [os.path.join('..', 'README.md')]},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.0.3',
    description='Simulation of rooted phylogenetic trees under a given Multitype Birth–Death model.',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://gitlab.pasteur.fr/phylo/treesimulator',
    keywords=['phylogenetics', 'tree generator', 'multitype birth-death model'],
    install_requires=['ete3', 'numpy'],
    entry_points={
            'console_scripts': [
                'generate_bd = treesimulator.simulate_forest_bd:main',
                'generate_bdei = treesimulator.simulate_forest_bdei:main',
                'generate_bdss = treesimulator.simulate_forest_bdss:main',
                'generate_mtbd = treesimulator.simulate_forest:main',
            ]
    },
)
