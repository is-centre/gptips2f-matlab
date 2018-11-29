# GPTIPS2F Toolbox for MATLAB #

GPTIPS2F is the evolution of the second version of the MATLAB toolbox developed by Dr. Dominic Searson. It is maintained by Dr. Aleksei Tepljakov, TalTech University, https://atdesign.ee/.

### New features ###

* Preset Random Constants (PRCs): a subset of Ephemeral Random Constants (ERCs) the difference being that PRCs are randomly chosen from a predefined set.
* Automatically Defined Functions (ADFs): basically templates that are seeded into the initial population and can also arise naturally during mutation.

### Installation ###

Run the code `addpath(genpath(gptips2f))` and save the resulting path afterwards.

### A brief on how to use ADFs ###

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
