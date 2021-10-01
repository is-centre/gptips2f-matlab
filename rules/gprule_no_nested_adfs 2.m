function [gp_out, fit] = gprule_no_nested_adfs(gp, expr, params)
%GPRULE_NO_NESTED_ADFS Individuals with nested ADFs are discarded
%   Individuals having the trait of being expressed through nested ADFs
%   should not have any chances of survival.

gp_out = gp;

% Params are not used at the moment
params; %#ok

% In expressions, ADFs are always encoded as certain known symbols
adfs = gp.nodes.adf.seed_str;

adf_found = false; % Found ADF
open_count = 0; % Open parentheses counter
k = 1; % Character pointer
fit = 1.0; % By default, the gene is fit for survival

% Parser body
% Here we assume that some rules for expression strings hold,
% e.g., that they are always complete [TODO: potential bug]
while k <= length(expr)
    
    if ~isempty(strfind(adfs, expr(k))) %#ok: backwards compatibility
        if adf_found % Nesting detected
            fit = 0.0;
            break;
        else
            adf_found = true;
            k = k + 1;
        end
    end
    
    if strcmpi(expr(k), '(') && adf_found
        open_count = open_count + 1;
    end
    
    if strcmpi(expr(k), ')') && adf_found
        open_count = open_count - 1;
        if open_count == 0
            adf_found = false;
        end
    end
    
    k = k + 1; % Increment character counter

end

end

