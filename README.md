# GPTIPS2F Toolbox for MATLAB #

GPTIPS2F is the evolution of the second iteration of the MATLAB toolbox developed by Dominic Searson.

### New features ###

* Preset Random Constants (PRCs): a subset of Ephemeral Random Constants (ERCs) the difference being that PRCs are randomly chosen from a predefined set.
* Automatically Defined Functions (ADFs): basically templates that are seeded into the initial population and can also arise naturally during mutation.
* Function nodes that allow to evolve solutions of (fractional-order) differential equations.

### Installation ###

* Run the code 
```
#!matlab

addpath(genpath(gptips2f))
```
and save the resulting path afterwards.