function [gp_out, fit] = gprule_no_nested_adfs(gp, expr, params)
%GPRULE_NO_NESTED_ADFS Individuals with nested ADFs are discarded
%   Individuals having the trait of being expressed through nested ADFs
%   with a max_nesting_level greater than allowed should be discarded.
%
% Example configuration:
%
% % Max nesting level config
% op = struct;
% op.max_nesting_level = 2;
%
% % General settings
% gp.evolution.rules.use = true;
% gp.evolution.rules.sets = {{@gprule_no_nested_adfs, op}};
% gp.evolution.rules.attempts = 500;

gp_out = gp;

% Maximum nesting level
max_nesting_level = 1;
if isfield(params, 'max_nesting_level')
    max_nesting_level = params.max_nesting_level;
end

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

%     % DEBUG
%     disp('Have the following expression:');
%     disp(expr);
%     disp('Have the following ADFs:');
%     disp(gp.nodes.adf.seed_str);
    
    if ~isempty(strfind(adfs, expr(k))) %#ok: backwards compatibility
        if adf_found && open_count >= max_nesting_level % Nesting detected
            fit = 0.0;
%             % DEBUG
%             fprintf(2, 'Found %d level nesting! Discarding invidivual.\n', open_count);
% 
%             disp('Have the following expression:');
%             disp(expr);
%             disp('Have the following ADFs:');
%             disp(gp.nodes.adf.seed_str);
% 
%             pause;
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

