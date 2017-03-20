function treestr = treegen(gp,fixedDepth)
%TREEGEN Generate a new encoded GP tree expression.
%
%   TREESTR = TREEGEN(GP) generates an encoded symbolic expression TREESTR
%   using the parameters stored in the structure GP.
%
%   TREESTR = TREEGEN(GP,DEPTH) generates an expression TREESTR of depth
%   DEPTH but otherwise using the parameters stored in the GP structure.
%
%   Trees are nominally built to a depth specified in the field
%   GP.TREEDEF.MAX_DEPTH according to the method specified in the field
%   GP.TREEDEF.BUILD_METHOD. However, if the optional additional input
%   argument DEPTH is specified (which must be an integer > 0) then a full
%   tree of this depth will be built irrespective of the settings in the GP
%   structure.
%
%   Remarks:
%
%   GPTIPS uses a encoded string representation of expressions which is not
%   very human readable, but is fairly compact and makes events like
%   crossover and mutation easier to handle. An example of such an
%   expression is a(b(x1,a(x4,x3)),c(x2,x1)). However, before an expression
%   is evaluated it is converted using TREE2EVALSTR to produce an evaluable
%   math expression, e.g. times(minus(x1,times(x4,x3)),plus(x2,x1)). In
%   GPTIPS, a tree of depth 1 contains a single node.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also EXTRACT, TREE2EVALSTR, GPREFORMAT

%Check input arguments for depth overide parameter
if nargin < 2
    fixedDepth = 0;
end

%extract parameters from gp structure
maxDepth = gp.treedef.max_depth;          %will use user max depth unless overridden by input argument
buildMethod = gp.treedef.build_method;    %1=full 2=grow 3=ramped 1/2 and 1/2
p_ERC = gp.nodes.const.p_ERC;             %terminal config. [0=no constants, 0.5=half constants half inputs 1=no inputs]
p_int = gp.nodes.const.p_int;             %probability that an ERC will be an integer
num_inp = gp.nodes.inputs.num_inp;
range = gp.nodes.const.range;
rangesize = gp.nodes.const.rangesize;
fi = gp.nodes.const.format;
ERCgenerated = false;
PRCgenerated = false;

afid_argt0  = gp.nodes.functions.afid_argt0; %functions with arity>0
afid_areq0  = gp.nodes.functions.afid_areq0; %functions with arity=0;
arity_argt0 = gp.nodes.functions.arity_argt0; %arities of functions
fun_lengthargt0 = gp.nodes.functions.fun_lengthargt0;
fun_lengthareq0 = gp.nodes.functions.fun_lengthareq0;

%if a fixed depth tree is required use 'full' build method
if fixedDepth
    maxDepth = max(fixedDepth,1);
    buildMethod = 1;
end

%if using ramped 1/2 and 1/2 then randomly choose max_depth and
%build_method
if buildMethod == 3
    maxDepth = ceil(rand * gp.treedef.max_depth);
    buildMethod = ceil(rand * 2);  %set either to 'full' or 'grow' for duration of function
end

% Initial string structure (nodes/subtrees to be built are marked with the $
%token). If ADFs are enabled, generate one of the ADFs from the list with
%given probability.
if gp.nodes.adf.use && rand <= gp.nodes.adf.p_gen
    % Choose an ADF to generate from the seed list
    treestr = ...
        ['(' ...
        gp.nodes.adf.use_seeds{ceil(rand*length(gp.nodes.adf.use_seeds))} ...
        ')'];
    
    % Make sure constants are actually generated at the end of tree creation
    if ~isempty(strfind(treestr,'?')), ERCgenerated = true; end
    if ~isempty(strfind(treestr,'#')), PRCgenerated = true; end
    
else
    treestr = '($)';
end

%recurse through branches
while true
    
    %find leftmost node token $ and store string position
    nodeTokenInd = strfind(treestr,'$');
    
    %breakout if no more nodes to fill
    if isempty(nodeTokenInd)
        break
    end
    
    %get next node position and process it
    nodeToken = nodeTokenInd(1);
    
    %replace this node token with 'active' node token, @
    treestr(nodeToken) = '@';
    
    %count brackets from beginning of string to @ to work out depth of @
    left_seg = treestr(1:nodeToken);
    numOpen = numel(strfind(left_seg,'('));
    numClosed = numel(strfind(left_seg,')'));
    depth = numOpen - numClosed;
    
    %choose type of node to insert based on depth and building method. If
    %root node then pick a non-terminal node (unless otherwise specified).
    if depth == 1 || (gp.nodes.adf.use && depth == 2)
        nodeType = 1; % 1=internal 2=terminal 3=either
        if maxDepth == 1 %but if max_depth is 1 must always pick a terminal
            nodeType = 2;
        end
        
        %if less than max_depth and method is 'full' then only pick a
        %function with arity>0
    elseif buildMethod == 1 && depth < maxDepth
        nodeType = 1;
        % TODO: temporary solution, must find out whether there is conflict
        % between this condition and the fulfilment of maximum tree depth
    elseif depth >= maxDepth % if depth is max_depth then just pick terminals
        nodeType = 2;
    else %pick either with equal prob.
        nodeType = ceil(rand * 2);
    end
    
    if nodeType == 1 %i.e a function with arity>0
        funChoice = ceil(rand * fun_lengthargt0);
        funName = afid_argt0(funChoice);
        numFunArgs = arity_argt0(funChoice);
        
        funRepStr = '';
        
        % For ADFs, check if argument enforcing is turned on
        if gp.nodes.adf.use && gp.nodes.adf.arg_force
            
            % Check whether this function is an ADF
            adfind = strfind(gp.nodes.adf.seed_str, funName);
            if ~isempty(adfind)
                funRepStr = gp.nodes.adf.use_seeds{adfind(1)};
                % Make sure constants are actually generated at the end of tree creation
                if ~isempty(strfind(funRepStr,'?')), ERCgenerated = true; end
                if ~isempty(strfind(funRepStr,'#')), PRCgenerated = true; end
            end
            
        end
        
        % Else just create a pattern from the $ terminals
        if isempty(funRepStr)
            %create appropriate replacement string e.g. a($,$) for 2 argument
            %function
            funRepStr = [funName '($'];
            if numFunArgs > 1
                for j=2:numFunArgs
                    funRepStr = [funRepStr ',$'];
                end
                funRepStr = [funRepStr ')'];
            else % for single argument functions
                funRepStr = [funName '($)'];
            end
        end
        
        %replace active node token @ with replacement string
        treestr = strrep(treestr,'@',funRepStr);
        
    elseif nodeType == 2 %pick a terminal (or an arity 0 function, if active)
        
        %choose a function or input with equal probability
        if fun_lengthareq0 && (rand >= 0.5 || (num_inp==0 && p_ERC==0))
            
            funChoice = ceil(rand * fun_lengthareq0);
            funName = afid_areq0(funChoice);
            
            %create appropriate replacement string for 0 arity function
            funRepStr = [funName '()'];
            
            %now replace active node token @ with replacement string
            treestr = strrep(treestr,'@',funRepStr);
            
        else %choose an input (if it exists) or
            %check if ERCs are switched on, if so pick a ERC
            %node/arithmetic tree token with p_ERC prob.
            term = false;
            
            if rand >= p_ERC
                term = true;
            end
            
            if term %then use an ordinary terminal not a constant
                inpChoice = ceil(rand * num_inp);
                treestr = strrep(treestr,'@',['x' sprintf('%d',inpChoice)]);
                
            else %use an ERC *or* PRC token (? or # character, respectively)
                % Which token to use?
                if rand <= gp.nodes.pconst.p_PRC
                    tok = '#';
                    PRCgenerated = true;
                else
                    tok = '?';
                    ERCgenerated = true;
                end
                treestr = strrep(treestr,'@',tok);
            end
        end
        
    end
    
end

%constant processing
if ERCgenerated
    
    %find ERC locations
    constInd = strfind(treestr,'?');
    numConstants = numel(constInd);
    
    %generate reqd. number of constants in the allowed range:
    const = rand(1,numConstants) * rangesize + range(1);
    
    %process integer ERCs if enabled
    if p_int
        ints = rand(1,numConstants) <= p_int;
        const(ints) = round(const(ints));
    else
        ints = zeros(1,numConstants);
    end
    
    %loop through ERC locations and replace with constant values
    for k=1:numConstants;
        
        if ints(k)
            arg = ['[' int2str(const(k)) ']'];
        else
            arg = ['[' sprintf(fi,const(k)) ']'];
            
        end
        %replace token with const
        main_tree = extract(constInd(1),treestr);
        treestr = strrep(main_tree,'$',arg);
        constInd = strfind(treestr,'?');
    end
    
end

% Preset random constant processing
if PRCgenerated
    
    %find PRC locations
    constInd = strfind(treestr,'#');
    numConstants = numel(constInd);
    PRCs = gp.nodes.pconst.set;
    numPRCs = length(PRCs);
    
    for k=1:numConstants
        myconst = PRCs(ceil(rand*numPRCs));
        arg = ['[' sprintf(fi,myconst) ']'];
    end
    
    % Loop through PRC locations and replace with constant values
    % that are randomly picked from a predefined set
    for k=1:numConstants;
        % Replace token with const
        main_tree = extract(constInd(1),treestr);
        treestr = strrep(main_tree,'$',arg);
        constInd = strfind(treestr,'#');
    end
    
end

%strip off outside brackets
treestr = treestr(2:end-1);