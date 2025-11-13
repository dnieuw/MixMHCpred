from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="MixMHCpred",
    version="3.1",
    author="David Nieuwenhuijse",
    author_email="d.nieuwenhuijse@erasmusmc.nl",
    description="A pan-Allele predictor for MHC-I peptide binding",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dnieuw/MixMHCpred",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "logomaker",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "MixMHCpred=MixMHCpred.cli:main",
        ],
    },
    include_package_data=True,
    keywords="MHC peptide binding prediction bioinformatics",
)
