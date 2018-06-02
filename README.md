This is compilation of code used for simulation and for generating figures in the paper: **Echo statistics associated with discrete scatterers: A tutorial on physics-based methods**. The authors of this paper are [Timothy K. Stanton](https://www.whoi.edu/profile/tstanton/) ([@timstanton99](https://github.com/timstanton99)), [Wu-Jung Lee](https://leewujung.github.io/) ([@leewujung](https://github.com/leewujung)), and [Kyungmin Baik](mailto:kbaik@kriss.re.kr) ([@nupho27](https://github.com/nupho27)). @leewujung and @nupho27 jointly contribute to this repository. The paper is being re-submitted after the first round of review and a PDF copy will be uploaded here once accepted.

The code are written in Matlab and many of them involve using the parallel computing toolbox to speed up Monte Carlo simulation.

Use `git clone https://github.com/leewujung/echo_stat_tutorial` to clone this repo to your local machine.

Since code from another repo [`broadband-echo-stats`](https://github.com/leewujung/broadband-echo-stats) is used as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to plot Figure 23, you also need to run

`git submodule init` to initilize your local configulation file, and 

`git submodule update` to fetch check out the appropriate submodule commit.


Alternatively, you can do this when cloning the repo by using

`git clone --recurse-submodules https://github.com/leewujung/echo_stat_tutorial`

