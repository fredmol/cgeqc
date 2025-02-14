from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

__version__ = "1.0.0"

setup(
    name="cgeqc",
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    url="https://bitbucket.org/genomicepidemiology/cgeqc",
    license="Apache-2.0",
    install_requires=[
        "weasyprint>=63.0",
        "Jinja2>=3.1.4",
        "matplotlib>=3.9.2"
    ],
    author="Frederik Duus Møller",
    author_email="freddu@food.dtu.dk",
    description="CGE Quality Control Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={
        "console_scripts": [
            "cgeqc=bin.cgeqc:main"
        ]
    },
)

