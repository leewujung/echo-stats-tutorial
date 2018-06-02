### Brief introduction

This is a compilation of code used for generating figures in the paper: **Echo statistics associated with discrete scatterers: A tutorial on physics-based methods**. The authors of this paper are [Timothy K. Stanton](https://www.whoi.edu/profile/tstanton/) ([@timstanton99](https://github.com/timstanton99)), [Wu-Jung Lee](https://leewujung.github.io/) ([@leewujung](https://github.com/leewujung)), and [Kyungmin Baik](mailto:kbaik@kriss.re.kr) ([@nupho27](https://github.com/nupho27)). @leewujung and @nupho27 jointly contribute to this repository. The paper is re-submitted to the Journal of the Acoustical Society of America (JASA) after the first round of review and a PDF copy will be uploaded here once accepted.

The code are tested with Matlab R2017b and use functions in the Signal Processing Toolbox ver7.5 and Parallel Computing Toolbox ver 6.11.

### To use the code

You can get a copy of the code in the following ways:
  1. Download the whole repo as a zip file using the "Clone or download" button on the top-right corner of this page
  - Use `git clone --recurse-submodules https://github.com/leewujung/echo_stat_tutorial` to clone this repo to your local machine.

  Note the `--recurse-submodules` option because code from another repo [`broadband-echo-stats`](https://github.com/leewujung/broadband-echo-stats) is used as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to plot Figure 23.

  If you forgot to use this option, run `git submodule init` to initialize your local configuration file, and `git submodule update` to check out the appropriate submodule commit.
