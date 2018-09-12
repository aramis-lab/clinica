from os.path import dirname, join
from setuptools import setup, find_packages
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

with open(join(dirname(__file__), 'clinica/VERSION'), 'rb') as f:
    version = f.read().decode('ascii').strip()


setup(
    name='clinica',
    version=version,
    url='http://clinica.run',
    description='Software platform for clinical neuroimaging studies',
    long_description=open('README.md').read(),
    author='ARAMIS Lab',
    maintainer='Clinica developers',
    maintainer_email='clinica-user@inria.fr',
    license='MIT license',
    packages=find_packages(exclude=('tests', 'tests.*')),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['clinica = clinica.cmdline:execute']
    },
    classifiers=[
        'Framework :: Clinica',
        'Development Status :: 0.1.0 - Dev',
        'Environment :: Console',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
