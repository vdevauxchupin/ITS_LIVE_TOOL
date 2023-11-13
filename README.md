# ITS_LIVE_TOOL

```{note}
add figures, images of processed data, data access gifs etc.
```

# About

`ITS_LIVE_TOOL` is a package designed to aid users working with the
[Inter-mission Time Series of Land Ice Velocity and Elevation](link)
(ITS_LIVE) dataset. The package provides functions for accessing data as
well as various methods to process ITS_LIVE observations. This notebook
will demonstrate various elements of the package and walk through the
steps of a typical workflow using ITS_LIVE_TOOL.

# Installation 

eventually, we hope to have a pip install. for now, install via:

`pip install pip+https://github.com/vdevauxchupin/ITS_LIVE_TOOL`

# Overview

## Data access + organization

We implement 3 classes of objects to store ITS_LIVE velocity data with the goal of making it easy and intuitive to keep track of the data you're working with and scale analysis. 

The illustration below provides a high-level overview of the main object classes within this package: 
[ITS_LIVE_TOOL objects](https://github.com/e-marshall/ITS_LIVE_TOOL/blob/main/figs/Image_20231113_113529_391.jpeg?raw=true)

### Interactive widget

### [`Glacier`](https://e-marshall.github.io/ITS_LIVE_TOOL/setup.html#glacier), [`Glacier_Point`](https://e-marshall.github.io/ITS_LIVE_TOOL/setup.html#glacier_point), and `Glacier_Centerline` objects

Depending on your purposes, create an object of the
[`Glacier`](https://e-marshall.github.io/ITS_LIVE_TOOL/setup.html#glacier),
[`Glacier_Point`](https://e-marshall.github.io/ITS_LIVE_TOOL/setup.html#glacier_point)
or `Glacier_Centerline` classes based on your selection from the
interactive map. These objects contain both raster and vector data and
are built on the Open Global Glacier Model (OGGM) and data access
scripts developed by the developers of the ITS_LIVE dataset. Check out
the docs to learn more!

## Data Processing 

We demonstrate and make available two processing routines. Be sure to
check out the accompanying [book]() and consider if either of these are
appropriate for your data and use case. Note that these methods are in
active development and thus should be considered in *beta* phase. Please
perform your own data inspection and due diligence if implementing these
methods.

### Inversion

Description – link to full description and examples

### Gaussian Process

Description – link to full description and examples


## How to use

## Citing

## Contact

## Acknowledgements 

## References

## Contributing