function [gp_out, fit] = gprule_no_nested_adfs(gp, expr)
%GPRULE_NO_NESTED_ADFS Summary of this function goes here
%   Individuals having the trait of being expressed through nested ADFs
%   should not have high chances of survival.
gp_out = gp;
fit = 1.0;

end

