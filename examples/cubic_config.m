function gp=cubic_config(gp)
%CUBIC_CONFIG Config file for multigene regression on a simple cubic polynomial.
%
%   GP = CUBIC_CONFIG(GP) generates a parameter structure GP that specifies
%   the GPTIPS run settings for multigene regression on data generated by
%   the cubic polynomial:
%
%         3         2
%   3.4 x1  + 2.9 x1  + 6.2 x1 + 0.75
%
%   Example:
% 
%   GP = RUNGP(@cubic_config) uses this configuration file to perform
%   symbolic regression on the cubic polynonial data and returns the
%   results in the data struct GP.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also SALUSTOWICZ1D_CONFIG, UBALL_CONFIG, RIPPLE_CONFIG,
%   REGRESSMULTI_FITFUN, RUNGP

%generate cubic polynomial data for train and test sets
%and use GPTIPS defaults for everything else
xtr = (-10: 0.2 : 10)';
ytr = 3.4 * xtr.^3 + 2.9*xtr.^2 + 6.2*xtr + 0.75;
gp.userdata.xtrain = xtr;
gp.userdata.ytrain = ytr;

xte = (-9.5 : 0.333 : 9.5)';
yte = 3.4 * xte.^3 + 2.9*xte.^2 + 6.2*xte + 0.75;
gp.userdata.xtest = xte;
gp.userdata.ytest = yte;
gp.userdata.name = 'Cubic poly';

%termination criterion
gp.fitness.terminate = true;
gp.fitness.terminate_value = 1e-6;

%function nodes
gp.nodes.functions.name = {'times','minus','plus','add3','mult3'};

