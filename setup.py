from setuptools import find_packages, setup

setup(
    author='Jana Hoffmann', 
    author_email='epicore_jana@family-hoffmann.de', 
    python_requires='>=3.12', 
    description='Compute core epitopes from multiple overlapping peptides.',
    license='MIT license',
    name='epicore', 
    url='https://github.com/AG-Walz/epicore',
    entry_points={
        'console_scripts': ['epicore=epicore_utils.epicore_main:main']
    },
    install_requires=[
        'biopython>=1.80',
        'click>=8.1',
        'matplotlib>=3.4',
        'numpy>=2',
        'pandas>=2',
        'pyyaml>=6.0.2'
    ],
    packages=find_packages(),
    classifiers = [
        'Programming Language :: Python :: 3', 
        'Programming Language :: Python :: 3.12'
    ]
)