# GPTIPS2F Toolbox for MATLAB #

GPTIPS2F is the evolution of the second version of the MATLAB toolbox developed by Dr. Dominic Searson.

Since 2017, a fork of the toolbox (‘2F’) is maintained by Dr. Aleksei Tepljakov, TalTech University, https://taltech.ee/

**Note that the version number of GPTIPS2F has been reset to 1.0 since Dec 6, 2018.**

## New features ##

* Preset Random Constants (PRCs): a subset of Ephemeral Random Constants (ERCs) the difference being that PRCs are randomly chosen from a predefined set.
* Automatically Defined Functions (ADFs): basically templates that are seeded into the initial population and can also arise naturally during mutation.
* Evolutionary rules: define rules for discarding individuals with certain undesirable traits or significantly decrease their chances of staying in the population (alpha testing feature).

## Installation ##

Run the code `addpath(genpath(gptips2f))` and save the resulting path afterwards.

## Documentation ##
Please consider [the original docs](https://sites.google.com/site/gptips4matlab/). Updated docs will be published in due time. Meanwhile, read on to learn how to use the ADF feature.

## A brief on how to use ADFs ##

Contrary to the original meaning of the term, ADFs---automatically defined functions---are defined manually by the user prior to the regression procedure. The purpose of ADFs at the moment is to enforce some structure into the otherwise unstructured modeling problem. Think of them as intelligent design elements---those that would be most useful in a certain context.

Here's how to set up ADFs in MATLAB code in the gp config file:
```
% ADF expressions. Can contain any custom function also.
gp.nodes.adf.expr = {'$1+?1*$2', '$1*(#1+$2)'};

% For terminals (the arguments of the ADF functions) use
%     $1, $2, ... to define normal terminals
%     ?1, ?2, ... to mark a terminal where ERC should go
%     #1, #2, ... to mark a terminal where PRC should go
% Don't forget to enumerate these correctly! Counters are unique to each ADF
% (see example above)
   
% Assign names to ADFs
gp.nodes.adf.name = {'simple_rnd_sum', 'simple_rnd_prod'};

% Predefined random constants (PRC) settings

% If a PRC is generated (see above), only the entries in this set are used
gp.nodes.pconst.set = [0.1 0.2 0.3 0.4 0.5]; 

% General probability of generating a PRC (does not affect "#nn" terminals)
gp.nodes.pconst.p_PRC = 0.25;

% Probability of generating an ADF
gp.nodes.adf.p_gen = 0.25;

% Enforce structure on mutation
gp.nodes.adf.arg_force = true;
```

## A brief on how to use Rules ##

Evolutionary rules allow to either discard or otherwise decrease the chances of survival of individuals having certain traits. These traits must be extracted from symbolic gene expressions and tested in user-defined rule functions. An example of how to structure the latter can be found in the `rules` folder under the name `gprule_template.m`.

There is also an implemented rule called `gprule_no_nested_adfs.m` available in the same folder. This rule will allow to remove from the population any individual who has genes where nested ADFs appear. This is useful for a certain class of regression problems. In the config function, this rule is set up as follows:
```
% Evolutionary rules
gp.evolution.rules.use = true;
gp.evolution.rules.sets = {{@gprule_no_nested_adfs, []}};
```

**Note that this is an alpha feature.**