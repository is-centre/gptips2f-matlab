function [symbolicResult, parallelResult, statsResult] = gptoolboxcheck
%GPTOOLBOXCHECK Checks if certain toolboxes are installed and licensed.
%
%   Checks for the Stats, Parallel Computing and Symbolic Math toolboxes.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPINIT

symbolicResult = @() (ismember('Symbolic Math Toolbox',gp_toolboxlist()) && license('test','symbolic_toolbox'));
parallelResult = @() (ismember('Parallel Computing Toolbox',gp_toolboxlist()) && license('test','distrib_computing_toolbox'));
statsResult = @() (( ismember('Statistics Toolbox',gp_toolboxlist()) || ismember('Statistics and Machine Learning Toolbox',gp_toolboxlist()) ) ...
    && license('test','statistics_toolbox')); 

