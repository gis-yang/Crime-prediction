# ST-Cokriging ArcGIS extension for crime prediction


## Introduction
**1. Geostatistics and Cokriging**

In geostatistics, Cokriging is a multivariate variant of Kriging technique and makes spatial predictions for a sparsely sampled variable (the primary variable) of interest, with help of one or more well-sampled ancillary variables (the secondary co-variables). Cokriging method usually results in more accurate predictions of the target primary variable than Kriging. This is because Cokriging method exploits cross-correlations between the primary variable and the secondary co-variables in addition to the spatial autocorrelation of the primary variable.
In conventional Cokriging method, both the primary variable and co-variables are in spatial domain, and time dimension is not taken into consideration. By extending it from sole spatial domain to the spatio-temporal domain, this algorithm formulated a ST-Cokriging method that is capable of taking advantage of both spatial and temporal correlation within and between primary and co-variables to produce temporally frequent predictions for the primary variable at a high spatial resolution as the co-variables.

<img src="/images/fg1.png">

Figure 1. Data processing flowchart of ST-Cokriging method 

**2.	ST-Cokriging work flow**

In ST-cokriging formulation, we assume that the primary variable of interest is coarse spatial resolution images that are sampled at a high temporal frequency (high temporal resolution), and the secondary variable (co-variable) are ﬁne spatial resolution images that are sparsely sampled over time (low temporal resolution), as shown in Figure 1. Without loss of generality, we only consider the case with only one co-variable observed at multiple time points in the mathematical formulation of ST-Cokriging method. The extension to two or more co-variables observed at multiple time points is straightforward.

**3.	Extension implementation**
In the Cokriging linear system, the covariance matrix C is of very large size. Since the co-variable is observed at high spatial resolution, it is very likely that M_j’s are very large. Similarly, even with smaller N_i’s for the primary variable observed at coarse spatial resolution, ∑_(i=1)^T▒N_i  can still be large, since the primary variable is observed at high temporal frequency. Therefore, the matrix C can be of high dimension. Solving such a high-dimensional linear system can be computationally infeasible. One popular method to alleviate this diﬃculty is to force small numbers in the matrix (vector) to be zero, known as thresholding or tapering. Meanwhile, we take advantage of the feature of regularly gridded data in the application presented in Section 3, which facilitates eﬃcient parallel computing of the Cokriging predictor and variance.

## Source and Environmental 

1.	Source code download
	- Download the toolbox in ArcGIS from here.
	- Download the script for spatio-temporal semi-variogram here.
	- Download the script for ST-Cokriging here.
2.	Environmental setup
	- Enable extensions in ArcMap
	<img src="/images/fg2.png">
	- Setup the Geo-processing option
	<img src="/images/fg3.png">
	- Open Arctoolbox window, right click and add a toolbox
	<img src="/images/fg4.png">
	- Navigate to the toolbox just downloaded and select
	<img src="/images/fg5.png">
	


A Github Pages template for academic websites. This was forked (then detached) by [Stuart Geiger](https://github.com/staeiou) from the [Minimal Mistakes Jekyll Theme](https://mmistakes.github.io/minimal-mistakes/), which is © 2016 Michael Rose and released under the MIT License. See LICENSE.md.

I think I've got things running smoothly and fixed some major bugs, but feel free to file issues or make pull requests if you want to improve the generic template / theme.

### Note: if you are using this repo and now get a notification about a security vulnerability, delete the Gemfile.lock file. 

# Instructions

1. Register a GitHub account if you don't have one and confirm your e-mail (required!)
1. Fork [this repository](https://github.com/academicpages/academicpages.github.io) by clicking the "fork" button in the top right. 
1. Go to the repository's settings (rightmost item in the tabs that start with "Code", should be below "Unwatch"). Rename the repository "[your GitHub username].github.io", which will also be your website's URL.
1. Set site-wide configuration and create content & metadata (see below -- also see [this set of diffs](http://archive.is/3TPas) showing what files were changed to set up [an example site](https://getorg-testacct.github.io) for a user with the username "getorg-testacct")
1. Upload any files (like PDFs, .zip files, etc.) to the files/ directory. They will appear at https://[your GitHub username].github.io/files/example.pdf.  
1. Check status by going to the repository settings, in the "GitHub pages" section
1. (Optional) Use the Jupyter notebooks or python scripts in the `markdown_generator` folder to generate markdown files for publications and talks from a TSV file.

See more info at https://academicpages.github.io/

## To run locally (not on GitHub Pages, to serve on your own computer)

1. Clone the repository and made updates as detailed above
1. Make sure you have ruby-dev, bundler, and nodejs installed: `sudo apt install ruby-dev ruby-bundler nodejs`
1. Run `bundle clean` to clean up the directory (no need to run `--force`)
1. Run `bundle install` to install ruby dependencies. If you get errors, delete Gemfile.lock and try again.
1. Run `bundle exec jekyll liveserve` to generate the HTML and serve it from `localhost:4000` the local server will automatically rebuild and refresh the pages on change.

# Changelog -- bugfixes and enhancements

There is one logistical issue with a ready-to-fork template theme like academic pages that makes it a little tricky to get bug fixes and updates to the core theme. If you fork this repository, customize it, then pull again, you'll probably get merge conflicts. If you want to save your various .yml configuration files and markdown files, you can delete the repository and fork it again. Or you can manually patch. 

To support this, all changes to the underlying code appear as a closed issue with the tag 'code change' -- get the list [here](https://github.com/academicpages/academicpages.github.io/issues?q=is%3Aclosed%20is%3Aissue%20label%3A%22code%20change%22%20). Each issue thread includes a comment linking to the single commit or a diff across multiple commits, so those with forked repositories can easily identify what they need to patch.
