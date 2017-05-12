from os.path import dirname, join
from setuptools import setup, find_packages

with open(join(dirname(__file__), 'clinica/VERSION'), 'rb') as f:
    version = f.read().decode('ascii').strip()

setup(
    name='Clinica',
    version=version,
    url='www.aramislab.fr',
    description='',
    long_description=open('README.md').read(),
    author='Aramis Lab',
    maintainer='Michael Bacci',
    maintainer_email='michael.bacci@inria.fr',
    license='TODO',
    packages=find_packages(exclude=('tests', 'tests.*')),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['clinica = clinica.cmdline:execute']
    },
    classifiers=[
        'Framework :: Clinica',
        'Development Status :: 0.1 - Dev',
        'Environment :: Console',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    install_requires=[
        'nibabel>=2.0.2',
        'nipype>=0.11.0,<0.13.0', # the t1-freesurfer pipeline does not run with nipype>=0.13.0
        'pybids>=0.1',
        'dipy>=0.6.0',
        'argcomplete>=1.4.1',
        'pandas>=0.18.1',
        'nose>=1.3.7'
    ],
    data_files=[('clinica', ['LICENSE.txt', 'clinica/VERSION'])],
)
