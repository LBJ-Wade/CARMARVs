This repo contains tools for running [CARMA
models](https://arxiv.org/abs/1802.09812) on radial velocity time series
containing up to one planetary signal and an arbitrary number of damped random
walk (real exponential correlation function) and damped harmonic oscillator
(complex exponential correlation function) stochastic components.

You will need to install some prerequisites in order to run the code here.  This
*should* be accomplished automatically by [Julia's](http://julialang.org) new
``Pkg`` mechanism; follow the example below.  It should be straightforward to
generalize to your particular data set and/or to script using Julia.

We have included an example data set an configuration file which is a RV
timeseries from [K2](https://keplerscience.arc.nasa.gov/) 39b, a system
containing one planet in a 4.6 day orbit and for which there is a claim of
asteroseismic power peaking at 340 muHz in [Van Eylen, et al.
(2016)](http://adsabs.harvard.edu/abs/2016AJ....152..143V).  We will fit the
planetary RV signal and model the asteroseismic power using a single oscillating
component.

First, we will need to download and install the prerequisites (note: in addition
to the Julia libraries that will be automatically downloaded, you will need the
Python libraries [corner](https://github.com/dfm/corner.py) and
[seaborn](https://seaborn.pydata.org/) for the visualizations here); begin by
initializing the environment.  Start Julia in the ``code`` directory (where you
find ``Project.toml``).  You can enter the ``Pkg`` REPL using the ``]`` key; do
so, and 'activate' the project:

```julia
(v1.0) pkg> activate .
```

If this is your first time running this project, then you will need to install
the prerequisites (you can skip this step once they are installed on your system
in subsequent runs):

```julia
(CARMARVs) pkg> instantiate
```

Now press backspace to exit the ``Pkg`` REPL back to the standard Julia REPL.
We import the code for running and examining CARMA models:

```julia
julia> using CARMARVs
```

The choices for various options running and post-processing of the code are
stored in a TOML formatted configuration file.  An example file with comments
can be found at ``K2-39b/runs/K2-39b.toml``; once you have defined the
parameters of the sampling, you can construct a posterior over CARMA models and
save the sampling using

```julia
julia> CARMARVs.run_carma_rvs("../K2-39b/runs/K2-39b.toml")
```

While it runs, this will occasionally checkpoint into the output directory; if
you (or your computer) kill the process before the termination criterion is
reached, then re-issuing the command will load the checkpoint file and resume
sampling.  If you want to really reset the sampling, just delete the checkpoint
file.

Once the sampling is done, the final state, including posteriors over the
Keplerian parameters, CARMA parameters, and calibration parameters for each
instrument (CSV file in the input directory) will be stored in the output
directory.  You can re-load them using

```julia
julia> posterior, samples, nest_state = CARMARVs.load_run("../K2-39b/runs/K2-39b.toml")
```

If you want to make the standard set of post-processing plots, including corner
plots of the Keplerian parameters, the CARMA parameters, a few draws from the
predicted RV state of the system given the data (including uncertainty),
inferred errorbar scalings, and a posterior over the inferred PSD, you can run

```julia
julia> CARMARVs.make_default_plots("../K2-39b/runs/K2-39b.toml")
```

Alternately, you can make your own calculations using the handy variable names
in the ``samples`` array.  For example, if you wanted to sigma-clip the
posterior over periods and plot the distribution (not that I advocate this!),
then you can run

```julia
julia> Periods = [s.P for s in samples]
julia> mu_P = mean(Periods)
julia> sigma_P = std(Periods)
julia> clipped = abs.((Periods .- mu_P) ./ sigma_P) .< 3
julia> sns.distplot(Periods[clipped])
```

See the ``RunCARMARVs.jl`` file for more examples of accessing the results of
your simulation.  Happy RV-ing!

If you would prefer to run with the parallel-tempered MCMC sampler from the
[Ensemble](https://github.com/farr/Ensemble.jl) package, see the example
notebook in the `code` directory.  This can be an order of magnitude more
efficient, but is not yet supported "out of the box."  Stay tuned!
