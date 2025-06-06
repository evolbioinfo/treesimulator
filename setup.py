from setuptools import setup, find_packages

setup(
    name='treesimulator',
    packages=find_packages(),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.2.18',
    description='Simulation of rooted phylogenetic trees under a given Multi-Type Birth–Death model (with or without Contact Tracing) (with or without a Skyline).',
    author='Anna Zhukova',
    author_email='anna.zhukova@pasteur.fr',
    url='https://github.com/evolbioinfo/treesimulator',
    keywords=['phylogenetics', 'tree generator', 'multitype birth-death model', 'contact tracing', 'skyline',
              'MTBD', 'MTBD-CT', 'MTBD-Skyline', 'MTBD-CT-Skyline'],
    install_requires=['six', 'ete3', 'numpy', 'scipy', 'sympy'],
    requires=['six', 'ete3', 'numpy', 'scipy', 'sympy'],
    entry_points={
            'console_scripts': [
                'generate_bd = treesimulator.simulate_forest_bd:main',
                'generate_bdei = treesimulator.simulate_forest_bdei:main',
                'generate_bdss = treesimulator.simulate_forest_bdss:main',
                'generate_bdeiss = treesimulator.simulate_forest_bdeiss:main',
                'generate_mtbd = treesimulator.simulate_forest:main',
            ]
    },
)
