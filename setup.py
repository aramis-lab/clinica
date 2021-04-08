import pathlib
from os.path import abspath, dirname, join

import pkg_resources
from setuptools import find_packages, setup

with pathlib.Path("requirements.txt").open() as requirements_txt:
    install_requires = [
        str(requirement)
        for requirement in pkg_resources.parse_requirements(requirements_txt)
    ]

with open(join(dirname(__file__), "clinica/VERSION"), "rb") as f:
    version = f.read().decode("ascii").strip()

this_directory = abspath(dirname(__file__))
with open(join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="clinica",
    version=version,
    url="http://www.clinica.run",
    description="Software platform for clinical neuroimaging studies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="ARAMIS Lab",
    maintainer="Clinica developers",
    maintainer_email="clinica-user@inria.fr",
    license="MIT license",
    packages=find_packages(exclude=("tests", "tests.*")),
    include_package_data=True,
    zip_safe=False,
    entry_points={"console_scripts": ["clinica = clinica.cmdline:execute"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Operating System :: OS Independent",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
    ],
    install_requires=install_requires,
    extras_require={"test": ["pytest", "coverage"]},
    python_requires=">=3.6,<3.8",
)
