
	evolecolpp - evolutionary ecology population genetic simulations 

# Background

This will be a collection of python modules for running evolutionary ecology simulations based on the [`fwdpp`](http://github.com/molpopgen/fwdpp) population genetics C++ template library.

# Modules

Currently there is only one example that implements the continuous snowdrift game modeled in [Doebeli, Hauert, and Killingback (2004, Science)](http://dx.doi.org/10.1126/science.1101456) and a metapopulation version of this game studied in [Wakano and Lehmann (2014, Journal of Theoretical Biology)](http://dx.doi.org/10.1016/j.jtbi.2014.02.036).

# Requirements

- python / numpy / matplotlib
- boost / boost-python
- fwdpp

# Installation

## Generic

After the dependencies, obtain the `evolecolpp` source and run `make`.

## Mac OS X

The easiest way to obtain the necessary packages is through [homebrew](http://github.com/Homebrew/homebrew).

## Linux

Installing the dependencies normally and install `fwdpp` from [source](http://github.com/molpopgen/fwdpp).

# Usage

An example simulation and plot can be run via `python cooperation_snowdrift.py` (compare to Figure 1A from Doebeli, Hauert, and Killingback (2004)):

![](https://raw.github.com/vancleve/evolecolpp/master/snowdrift_branching.png)

Or via `python cooperation_snowdrift_metapop.py` (similar to Figure 1A from Wakano and Lehmann (2014) though with different population size):

![](https://raw.github.com/vancleve/evolecolpp/master/snowdrift_branching_metapop.png)
