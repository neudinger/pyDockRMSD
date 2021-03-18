# DockRMSD

## Descriptions

[![GitHub license](https://img.shields.io/badge/license-EUPL-blue.svg)](https://raw.githubusercontent.com/herotc/hero-rotation/master/LICENSE) [![Build Github Status](https://github.com/neudinger/pyDockRMSD/workflows/Build%20pydockrmsd/badge.svg)](https://github.com/neudinger/pyDockRMSD/actions)

Docked Root-mean-square deviation of atomic positions [Paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0362-7)

![formula](https://render.githubusercontent.com/render/math?math={\mathrm{RMSD}=\sqrt{\frac{1}{N}\sum_{i=1}^N\delta_i^2}})

<!-- $$
\mathrm{RMSD}=\sqrt{\frac{1}{N}\sum_{i=1}^N\delta_i^2}
$$ -->

![formula](https://render.githubusercontent.com/render/math?math={\mathrm{RMSD}(\mathbf{v},\mathbf{w})=\sqrt{\frac{1}{n}\sum_{i=1}^n\|vi-w_i\|^2}=\sqrt{\frac{1}{n}\sum{i=1}^n((v{ix}-w{ix})^2+(v{iy}-w{iy})^2+(v{iz}-w{iz})^2})})

<!-- $$
\mathrm{RMSD}(\mathbf{v}, \mathbf{w})
= \sqrt{\frac{1}{n}\sum_{i=1}^n \|v_i - w_i\|^2}
= \sqrt{\frac{1}{n} \sum_{i=1}^n ((v_{ix} - w_{ix})^2 + (v_{iy} - w_{iy})^2 + (v_{iz} - w_{iz})^2})
$$ -->

## Paper extract

> Computer-aided drug design, in particular protein–ligand docking, has brought about the discovery of many biologically active drugs [1, 2]. In many protein–ligand docking programs, a flexible small molecule structure is docked in a rigid protein receptor structure in order to find the optimal binding conformation and affinity of the small molecule within the protein binding pocket. Since the ability of these programs to accurately assess binding affinity is dependent on their ability to find the optimal conformation of the ligand in the protein binding pocket, docking programs are often benchmarked by their ability to reproduce the native binding pose of a ligand from a protein–ligand complex crystal structure. A common metric used to evaluate distance between the predicted pose and the native pose, given a superposition of their protein receptor structures, is the root mean square deviation (RMSD) between their respective atoms

[DockRMSD PDF paper](https://zhanglab.ccmb.med.umich.edu/DockRMSD/DockRMSD.pdf)

## Usage

### Local Install

Linux only:

```bash
conda env create --name pydockrmsd --file condaenv/requirement.yml
conda activate pydockrmsd
conda env update --name pydockrmsd --file condaenv/ci-cd.yml
./scripts/install.sh
```

### Requirement

Build Requirement

```bash
pip install -r requirements.txt
```

Direct Download

```bash
pip install pydockrmsd # pypi source
```

## Example

- [crystal_bench.py](https://github.com/neudinger/pyDockRMSD/blob/main/examples/crystal_bench.py)
  - `rdkit, pandas` required
- [crystal_bench.ipynb](https://github.com/neudinger/pyDockRMSD/blob/main/examples/crystal_bench.ipynb)
[DockRMSD Website](https://zhanglab.ccmb.med.umich.edu/DockRMSD/) python Wrapper: docking pose distance calculation

Python atom mapping and RMSD calculation of symmetric molecules through graph isomorphism. ___Journal of Cheminformatics___, 11:40 (2019).

PyDockRMSD Written by Barre Kevin, DockRMSD Written by Eric Bell

## Os supported

- Linux
- Mac OS
- Window (soon)

Tools used:

- cython
- pdoc3
- pandoc
- cibuildwheel (cross compilation)

## Documentation

```bash
pdoc3 pydockrmsd --http localhost:8080
```

Tag:

- v: Version
- d: Documentation
- p: Publish
- t: Test
