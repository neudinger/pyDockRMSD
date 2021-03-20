import os
import pathlib
import re
from typing import List
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from Cython.Compiler import Options

Options.docstrings = True
Options.embed_pos_in_docstring = True
Options.embedsignature = True

here = pathlib.Path(__file__).parent.resolve()

version = os.getenv("version", default="0.0.0")
version = re.match(r".*?(\d+\.?\d+\.?\d+).*?",
                   version, re.S).groups()[0]

long_description = (here / 'README.md').read_text(encoding='utf-8')

print(f"__version__ = '{version}'",
      file=open(f"{here}/pydockrmsd/__version__.py", "w"))

extensions: List[Extension] = [
    Extension(
        name="pydockrmsd.dockrmsd",
        # Cannot be use on cross platform
        # extra_compile_args=["-static-libgcc", "--static", "-O3"],
        # extra_link_args=["-lm"],
        include_dirs=['pydockrmsd/DockRMSD_sources'],
        sources=["pydockrmsd/dockrmsd.pyx"],
    )
]

setup(
    name="pydockrmsd",
    version=version,
    author='Barre Kevin',
    packages=find_packages(),
    author_email='kevin.barre@epitech.eu',
    description='Python Dock RMSD calculation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # url='https://www.python.org/sigs/distutils-sig/',
    ext_modules=cythonize(extensions, annotate=False),
    python_requires=">=3.6",
    install_requires=[],
    classifiers=[
        # How mature is this project ? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        "Programming Language :: C",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Software Development",
        "Topic :: Scientific/Engineering",
        "Typing :: Typed",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    project_urls={  # Optional
        'Source': 'https://github.com/neudinger/pyDockRMSD',
    },
)
# os.remove(f"{here}/pydockrmsd/dockrmsd.c")
