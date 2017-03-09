function [fcn, newexpr, args, seed] = parseadf(expr, opt)
%PARSEADF Parse the ADF expression and return anonymous fcn, args, and seed
%
% Calling sequence:
% [FCN, NEWEXPR, ARGS, SEED] = PARSEADF(EXPR, OPT)
%
% where FCN contains the anonymous function for the direct evaluation of the
% expression, NEWEXPR contains the new expression string with generated
% arguments and ARGS contains a string containing these arguments encolsed
% in ( ). SEED contains the gene seed for generating the initial population.
% The input argument OPT is a string that contains either 's' or 'f'. In
% case of 's' the function outputs FCN as string (default), and in case of
% 'f' an anonymous function handle is returned.
%
% When designing your expressions, use any valid functions with the
% following symbols representing input arguments:
%
% $ - Free connection point
% # - Preset random constant (PRC) point that later turns into an ERC
% ? - An ephemeral random constant (ERC) point
%
% For example: 'times(#,plus($,?))'
%
% Note that the function must accept at least one argument. Depending on
% the context you can include constants into function calls. In the
% algorithm, the insides of the ADF are never revealed, so it is treated as
% just another function, with the big difference that you can specify
% locations to which terminal nodes (PRCs and ERCs) are attached with 100%
% probability during the initial creation of the population.
%
%   (c) Aleksei Tepljakov 2017
%
%   GPTIPS2F

% Check input
if ~isa(expr,'char')
    error('Expression must be a character array describing an ADF function.');
end

if nargin < 2
    opt = 's';
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

% Assign correct output
switch lower(opt)
    case 's'
        fcn = myfun;
    case 'f'
        fcn = eval(myfun);
    otherwise
        error('Unknown option given in second argument.');
end

end

