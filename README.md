[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1313729.svg)](https://doi.org/10.5281/zenodo.1313729)

## Introduction

This is a compilation of code used for generating figures in the following tutorial:

Stanton TK, Lee W-J, Baik K. (2018). **Echo statistics associated with discrete scatterers: A tutorial on physics-based methods**. The Journal of the Acoustical Society of America, 144(6): 3124–3171. [https://doi.org/10.1121/1.5052255](https://doi.org/10.1121/1.5052255).

**This tutorial is open-access without restriction directly from the JASA website**. Please email the authors if you have trouble accessing the complete tutorial.

[Wu-Jung Lee](https://leewujung.github.io/) ([@leewujung](https://github.com/leewujung)) and [Kyungmin Baik](mailto:kbaik@kriss.re.kr) ([@nupho27](https://github.com/nupho27)) jointly contributed to this repository. The code here are subject to future revisions while a frozen version of the code at the time of publication is available as [supplementary materials](https://asa.scitation.org/doi/suppl/10.1121/1.5052255) at the JASA website.


## To use the code

### Getting a copy
You can get a copy of the code in the following ways:
1. [Preferred] Use `git clone --recurse-submodules https://github.com/leewujung/echo_stat_tutorial` to clone this repo to your local machine.

    Note the `--recurse-submodules` option because code from another repo [`broadband-echo-stats`](https://github.com/leewujung/broadband-echo-stats) is used as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to plot Figure 23.

    If you forgot to use this option, run `git submodule init` to initialize your local configuration file, and `git submodule update` to check out the appropriate submodule commit.

2. Download the whole repo as a zip file using the "Clone or download" button on the top-right corner of this page.

    Note that the folder `broadband-echo-stats` will be empty after you unzip. This is because the folder was a git submodule linked to [a separate repo](https://github.com/leewujung/broadband-echo-stats). You can download the content of this folder [via this link](https://github.com/leewujung/broadband-echo-stats/archive/master.zip) and put everything under folder `broadband-echo-stats`. (As you can see, using the `git submodule` utility in 1. below is much more straightforward.)

### A few notes
1. The code are tested with Matlab R2017b and use functions in the Signal Processing Toolbox ver7.5 and Parallel Computing Toolbox ver 6.11.

2. Outputs from `figure_12.m` are used in many later code, so make sure to run this first. It will produce a folder `./figs/figure_12` which contains beampattern parameters and impulse response for a few different ka values that will be used for plotting other figures.

3. Many of these take _a while_ to run because the Monte Carlo simulation samples are generated on the fly by default (and we are taking a brute force approach here  without fancy sampling for simplicity). You can opt out of the simulation part by changing a flag (`mc_opt`) in the code if you have previously-calculated samples already. This applies to the following figures: 13, 15-18, 20, 21, and 23.

4. `figure_13` uses results from `figure_15`, so don't be alarmed if you get a warning asking you to run it first.

5. In the header of the code that generate figures 12-18, 20, 21, and 23, the note "3D distribution of scatterers" indicates that the scatterers are randomly and uniformly distributed within a thin hemispherical shell (half-space), and the note "2D distribution of scatterers" indicates that the scatterers are randomly and uniformly distributed within a thin arc of constant radius in the half-plane containing the MRA of the beam. See Sec. VII.A.2. in the paper for  detail.
