from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="archground",
    version="1.0.0",
    author="Juan Alday",
    description="Python package for ground-based planetary spectroscopy with archnemesis",
    long_description=long_description,
    long_description_content_type="text/markdown",  # important for Markdown rendering
    url="https://github.com/juanaldayparejo/archground-dist.git",
    project_urls={
        "Source": "https://github.com/juanaldayparejo/archground-dist.git",
    },
    #packages=["archnemesis"],
    packages=find_packages(), 
    install_requires=[
      'archnemesis',
      'numpy',
      'matplotlib',
      'numba>=0.57.0',
      'scipy',
      'joblib',
      'h5py',
      'astropy',
      'astroquery',
    ],
    extras_require={
        'docs': ['sphinx', 'sphinx_rtd_theme'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
