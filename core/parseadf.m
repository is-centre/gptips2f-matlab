function [fcn, seed, newexpr, args] = parseadf(expr)
%PARSEADF Parse the ADF expression and return an anonymous fcn and arg seed
%
% Calling sequence:
% [FCN, SEED, NEWEXPR, ARGS] = PARSEADF(EXPR)
%
% where NEWEXPR contains the new expression string with generated arguments
% and ARGS contains a string containing these arguments encolsed in ( ).
% 
%
% When designing your expressions, use any valid functions with the
% following symbols representing *all* input arguments:
%
% $ - Free connection point
% # - Preset random constant (PRC) point that later turns into an ERC
% ? - An ephemeral random constant (ERC) point
%
% For example: 'times(#,plus($,?))'
%
% Note that the use of constant numbers with the arguments is forbidden.
% In addition, the function must accept at least one argument.
%
%   (c) Aleksei Tepljakov 2017
%
%   GPTIPS 2F

% Check input
if ~isa(expr,'char')
    error('Expression must be a character array describing an ADF function.');
end

% Parse the arguments of the array
argname = 'x'; % Default argument is "x"
argind = 1;    % The index of the argument
seed = '';     % The seed is used when the initial population is built
exprind = 1;   % The point from which to take a piece of the old expr
newexpr = '';  % Updated expression with proper arguments inserted
for k=1:length(expr)
    switch expr(k)
        case {'$', '#', '?'}
            seed = [seed expr(k) ','];
            newexpr = [newexpr expr(exprind:(k-1)) argname num2str(argind)];
            exprind = k+1; argind = argind + 1;
        otherwise
            % Do nothing
    end
end

% Finalize the expression
if exprind <= length(expr)
    newexpr = [newexpr expr(exprind:end)];
end

% Are there any arguments?
if argind == 1
    error('The ADF must have at least one input argument.');
end

% Amend the seed
seed = ['(' seed(1:end-1) ')'];

% Create the anonymous function
args = '';
for k=1:(argind-1)
   args = [args argname num2str(k) ','];
end

% Amend the arguments and construct the expression
args = ['(' args(1:end-1) ')'];
myfun = ['@' args ' ' newexpr];

% Finally, return the anonymous function handle
fcn = eval(myfun);

end

