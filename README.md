
	evolecolpp - evolutionary ecology population genetic simulations 

# Background

This will be a collection of python modules for running evolutionary ecology simulations based on the [`fwdpp`](http://github.com/molpopgen/fwdpp) population genetics C++ template library.

# Modules

Currently there is only one example that implements the continuous snowdrift game modeled in [Doebeli, Hauert, and Killingback (2004, Science)](http://dx.doi.org/10.1126/science.1101456).

# Requirements

- python / numpy / matplotlib
- boost / boost-python
- fwdpp

# Installation

## Generic

After the dependencies, obtain the `evolecolpp` source and run `make`. An example simulation and plot can be run via `python cooperation_snowdrift.py`.

## Mac OS X

The easiest way to obtain the necessary packages is through [homebrew](http://github.com/Homebrew/homebrew).

## Linux

Installing the dependencies normally and install `fwdpp` from [source](http://github.com/molpopgen/fwdpp).

