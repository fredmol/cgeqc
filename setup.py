from setuptools import setup, find_packages
from pathlib import Path
import re

# Read version from version.py
with open('cgeqc/version.py', 'r') as f:
    version_file = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        __version__ = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string in cgeqc/version.py")

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="cgeqc",
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'cgeqc': [
            'templates/*.html',
            'assets/*.*',
        ],
    },
    url="https://bitbucket.org/genomicepidemiology/cgeqc",
    license="Apache-2.0",
    install_requires=[
        "weasyprint>=63.0",
        "Jinja2>=3.1.4",
        "matplotlib>=3.9.2",
        "scipy>=1.10.0",
        "numpy>=1.24.0"
    ],
    author="Frederik Duus Møller",
    author_email="freddu@food.dtu.dk",
    description="CGE Quality Control Tool for sequence data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    scripts=["bin/cgeqc"],
)