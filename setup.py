from os.path import dirname, join
from setuptools import setup, find_packages

with open(join(dirname(__file__), 'clinica/VERSION'), 'rb') as f:
    version = f.read().decode('ascii').strip()

setup(
    name='Clinica',
    version=version,
    url='http://xxxx',
    description='',
    long_description=open('README.md').read(),
    author='Aramis Equipe',
    maintainer='Michael Bacci',
    maintainer_email='michael.bacci@inria.fr',
    license='MIT',
    packages=find_packages(exclude=('tests', 'tests.*')),
    include_package_data=True,
    zip_safe=False,
    # entry_points={
    #     'console_scripts': ['clinica = clinica.cmdline:execute']
    # },
    classifiers=[
        'Framework :: Clinica',
        'Development Status :: 0.1 - Dev',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    install_requires=[
        'nipype>=0.11.0'
    ],
)
