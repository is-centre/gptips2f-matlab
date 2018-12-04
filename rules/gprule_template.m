function [gp_out, fit] = gprule_template(gp, expr, params)
%GPRULE_TEMPLATE Use this template to define a custom rule for fitness
% These rules are applied to genes, no complete expressions
%
% Input arguments: GP is the main structure, EXPR is the encoded gene
% expression, PARAMS are additional parameters needed for running the rule
%
% Output arguments:  GP_OUT is the GP structure (can be left unassigned if
% unmodified) and FIT is the fitness of this particular gene as determined
% by the rule logic.
%
% NB! FIT has a completely different meaning here compared to the overall
% fitness function. Here, FIT = 0.0 means the individual is unfit for
% survival, while a FIT = 1.0 means the opposite.

% Params are left unused in this case (no params)
params; %#ok

% The structure is here if it is needed;
% NB! Do not forget to pass it as output!
gp_out = gp;

% This is custom decision logic, but the idea is that
% the output of decision must be constrained to [0,1].

% As implemented in this particular example, FIT is always "true" or "pass"
% depending on the context, unless the input expression is empty.
%
% However, the logic of the rule may be quite complicated.

fit = 1.0;
if isempty(expr)
    fit = 0.0;
end

end

