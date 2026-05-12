# PMD - Planar Multibody Dynamics Open-Source Software

PMD is a Python library for dynamic analysis, which allows the construction of 
mechanical models and their numerical simulation.

The source code is available at [GitHub](https://github.com/CangiGia/PMD).

::::{grid} 3
:gutter: 2

:::{grid-item-card} 
:class-body: text-center
:link: getting_started/overview
:link-type: doc
```{image} /_static/getting_started.svg
:height: 50
```

### Getting Started 

New to PMD? Check out the getting started guides. They contain an
introduction to PMD main concepts and links to additional tutorials.
:::

:::{grid-item-card} 
:class-body: text-center
:link: user_guide/user_guide
:link-type: doc
```{image} /_static/user_guide.svg
:height: 50
```

### User Guide

The user guide provides in-depth information on the
key concepts of PMD with useful background information and explanation.
:::

:::{grid-item-card} 
:class-body: text-center
:link: api/pmd.core
:link-type: doc

### API Reference

Complete reference for all public classes, joints, forces, and the
simulation solver.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2
:caption: Getting Started

getting_started/overview
getting_started/install
getting_started/first_model
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: User Guide

user_guide/user_guide
user_guide/bodies_markers
user_guide/joints
user_guide/forces
user_guide/drivers_functions
user_guide/units
user_guide/solvers
user_guide/analyses
user_guide/gui
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: API Reference

api/pmd.core
api/pmd.gui
```
