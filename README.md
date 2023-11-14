# ITS_LIVE_TOOL

```{note}
add figures, images of processed data, data access gifs etc. to readme 
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

### Temporal baseline thresholds

This is not necessary in all situations but can be a useful step to familiarize yourself with the dataset with which you're working. See the notebook titled `temporalbaseline_filter.ipynb` for more discussion and examples of this. 

### Gaussian Process

Gaussian Process models are Bayesian, non-parametric methods that can be helpful for solving regression problems. Gaussian Process models require some understanding of the characteristics of your dataset. In essence, we create a model by specifying a prior distribution and covariance functions that we believe reasonably match our dataset. We then train the model on the dataset and use the model to make predictions over a specified prediction-space. Because this is a Bayesian approach, the model predictions also contain uncertainty quantification. The uncertainty quantification refers to the range of possible posterior functions that fit the observations. We use the mean of all posterior functions as the model predictions. So when we describe the model predictions, what we're really referring to is the mean of all posterior functions from the Gaussian Process model. 



### Inversion

Description – link to full description and examples


## How to use

For a detailed, code-based walk through of ITS_LIVE_TOOL functionality, check out the roadmap.ipynb notebook. This will demonstrate how to access ITS_LIVE data and use the convenience data structures provided in this package. It will also briefly describe and show sample outputs from the processing methods implemented in this package. For in-depth explanations of the processing methods, check out the individual notebooks where those steps will be described in much more detail.

## Citing

- idk if its appropriate/required here but I see some packages have a readme section for how to cite them if someone uses their work ? 

## Contact
add our emails

## Acknowledgements 
- funding sources, advisors, other collaborators, anyone/thing else? 

## References
### Data
- cite ITS_LIVE, RGI v7 outlines, OGGM centerlines
### Methods
- cite gaussian process, temporal inversion literature here?
### Software
- cite python packages used 


## Contributing
maybe add something like: We welcome community contributions to this work! Please don't hesitate to raise an issue, start a discussion our reach out to us over email

