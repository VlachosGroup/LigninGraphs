#
# setup.py
#
# Installation script to get setuptools to install nextorch into
# a Python environment.
#
import os
import sys
import setuptools

# Import the lengthy rich-text README as the package's long
# description:
root_dir = os.path.dirname(__file__)

with open(os.path.join(root_dir, "README.rst"), "r") as fh:
	long_description = fh.read()


setuptools.setup(
    name="ligning", 
    version="0.0.1",
    author="Vlachos Research Group",
    author_email="vlachos@udel.edu",
    description="Accelerated Lignin Structure Generation with Graph-based Multiscale Modeling ",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/VlachosGroup/ligningraphs",
    project_urls={
        "Documentation": "https://ligningraphs.readthedocs.io/en/latest/",
        "Source": "https://github.com/VlachosGroup/ligningraphs",
    },
    packages=setuptools.find_packages(),
    package_data={'':['*.xlsx']},
    include_package_data=True,
    python_requires=">=3.7",
    install_requires=[
        "matplotlib>=3.1.1",
        "numpy>=1.19.2",
        "scipy>=1.3.1",
        "pandas>=0.25.1",
        "openpyxl>=3.0.7",
        "pytest>=6.2.3",
        "networkx>=2.5",
        "pysmiles>=1.0.1",
        "rdkit-pypi>=2021.9.2.1"],
    classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Chemistry",
    ],
)