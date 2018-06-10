## Brief introduction

This is a compilation of code used for generating figures in the paper: **Echo statistics associated with discrete scatterers: A tutorial on physics-based methods**. The authors of this paper are [Timothy K. Stanton](https://www.whoi.edu/profile/tstanton/) ([@timstanton99](https://github.com/timstanton99)), [Wu-Jung Lee](https://leewujung.github.io/) ([@leewujung](https://github.com/leewujung)), and [Kyungmin Baik](mailto:kbaik@kriss.re.kr) ([@nupho27](https://github.com/nupho27)). @leewujung and @nupho27 jointly contribute to this repository. The paper is re-submitted to the Journal of the Acoustical Society of America (JASA) after the first round of review and a PDF copy will be uploaded here once accepted.

## To use the code

You can get a copy of the code in the following ways:
  1. Download the whole repo as a zip file using the "Clone or download" button on the top-right corner of this page
  2. Use `git clone --recurse-submodules https://github.com/leewujung/echo_stat_tutorial` to clone this repo to your local machine.

      Note the `--recurse-submodules` option because code from another repo [`broadband-echo-stats`](https://github.com/leewujung/broadband-echo-stats) is used as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to plot Figure 23.

      If you forgot to use this option, run `git submodule init` to initialize your local configuration file, and `git submodule update` to check out the appropriate submodule commit.

### A few things to note:
1. The code are tested with Matlab R2017b and use functions in the Signal Processing Toolbox ver7.5 and Parallel Computing Toolbox ver 6.11.

2. Outputs from `figure_12.m` are used in many later code, so make sure to run this first.

3. Many of these take _a while_ to run because the Monte Carlo simulation samples are generated on the fly by default (and we are taking a brute force approach here for simplicity). You can opt out the simulation part by changing a flag (`mc_opt`) in the code to plot pre-calculated samples. This applies to the following figures: 13, 15-18, 20-21, 23.

4. `figure_13` uses results from `figure_15`, so don't be alarmed if you get a warning asking you to run it first.


