function [gp_out, fit] = gprule_template(gp, expr)
%GPRULE_TEMPLATE Use this template to define a custom rule for fitness

% The structure is here if it is needed;
% NB! Do not forget to pass it as output!
gp_out = gp;

% This is custom decision logic, but the idea is that
% the output of decision must be constrained to [0,1].

% Here, FIT is always "true" or "pass" depending on the context.
% This means that the expression passed in EXPR is always evaluated with
% max fitness, unless it is empty.
%
% However, the logic of the rule may be quite complicated.

fit = 1.0;
if isempty(expr)
    fit = 0.0;
end

end

