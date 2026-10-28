from pathlib import Path

from setuptools import setup


ROOT = Path(__file__).parent

setup(
    name="MeHscr",
    version="2.1.0",
    description="Profile methylation heterogeneity from bisulfite-sequencing data",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Ya-Ting Sabrina Chang",
    author_email="ytchang.sabrina@gmail.com",
    url="https://github.com/britishcoffee/MeH",
    license="MIT",
    py_modules=["MeHscr"],
    python_requires=">=3.8",
    install_requires=[
        "joblib>=1.0",
        "numpy>=1.20",
        "pandas>=1.5",
        "pysam>=0.16",
    ],
    entry_points={"console_scripts": ["mehscr=MeHscr:main"]},
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
