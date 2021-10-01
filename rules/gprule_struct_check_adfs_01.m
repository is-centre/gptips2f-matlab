function [gp_out, fit] = gprule_struct_check_adfs(gp, expr, params)
%GPRULE_NO_NESTED_ADFS Individuals with nested ADFs are discarded
%   Individuals having the trait of being expressed through nested ADFs
%   should not have any chances of survival.

gp_out = gp;

% Params are not used at the moment
params; %#ok

% In expressions, ADFs are always encoded as certain known symbols
%adfs = gp.nodes.adf.seed_str;
%adf_found = false; % Found ADF
 % Character pointer
fit = 1.0; % By default, the gene is fit for survival
k=1;
%Reg expression
regExpr='(x1,x2,?';


% Arrays


% Here we assume that some rules for expression strings hold,
% e.g., that they are always complete [TODO: potential bug
   
if isempty(regexp(expr, regExpr))
        %disp(expr(k));
        %if adf_found % Nesting detected
   fit = 0.0;
   %disp('killed this');
   disp(expr);
else
    %disp(expr);
    
end
    
end