from setuptools import find_packages, setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

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
    install_requires=requirements,
    packages=find_packages(),
    classifiers = [
        'Programming Language :: Python :: 3', 
        'Programming Language :: Python :: 3.12'
    ]
)