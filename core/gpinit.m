function gp = gpinit(gp)
%GPINIT Initialises a run.
%
%   GP = GPINIT(GP)
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPCHECK, GPDEFAULTS, GPINITPARALLEL

%determine status of symbolic, parallel and stats toolboxes
[gp.info.toolbox.symbolic, gp.info.toolbox.parallel, gp.info.toolbox.stats] = gptoolboxcheck;

% Preprocess function nodes to account for ADFs
gp = procadfs(gp);

%process function nodes before run
gp = procfuncnodes(gp);

%throw an error if there are no inputs, p_ERC=0 and there are no arity zero
%functions active
if  gp.nodes.inputs.num_inp == 0 && gp.nodes.const.p_ERC == 0 && ...
        isempty(find(gp.nodes.functions.arity(logical(gp.nodes.functions.active)) == 0, 1))
    error('No terminals (inputs, constants or zero arity functions) have been defined for this run.');
end

%initialise some state and tracker variables
gp.state.count = 1;
gp.state.init_val = false; % AT: fixes bug where for some reason valbest is not generated on first generation
gp.state.best.fitness = [];
gp.state.best.individual = [];
gp.state.run_completed = false;
gp.state.current_individual = [];
gp.state.std_devfitness = [];
gp.state.terminate = false;
gp.state.force_compute_theta = false;
gp.fitness.returnvalues = cell(gp.runcontrol.pop_size,1);

%process mutation probabilities vector
gp.operators.mutation.cumsum_mutate_par = cumsum(gp.operators.mutation.mutate_par);

%init. history variables
gp.results.history.bestfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.meanfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.std_devfitness = zeros(gp.runcontrol.num_gen,1);
gp.results.history.about = 'Fitness on training data';
gp.results.history = orderfields(gp.results.history);

%best of run (on training data) fields
gp.results.best.fitness = [];
gp.results.best.individual = [];
gp.results.best.returnvalues = [];
gp.results.best.foundatgen = [];
gp.results.best.about = 'Best individual on training data';
gp.results.best = orderfields(gp.results.best);

%assign field holding fitnesses to gp structure
gp.fitness.values = zeros(gp.runcontrol.pop_size,1);
gp.fitness.complexity = zeros(gp.runcontrol.pop_size,1);

if strncmpi(func2str(gp.fitness.fitfun),'regressmulti',12)  && ~isfield(gp.userdata,'bootSampleSize')
    gp.userdata.bootSampleSize = size(gp.userdata.ytrain,1);   
end

%cache init
if gp.runcontrol.usecache
    gp = initcache(gp);
end

if ~gp.runcontrol.quiet
    fns = [];
    for i=1:length(gp.nodes.functions.active_name_UC);
        fns = [fns ' ' gp.nodes.functions.active_name_UC{i}];
    end
    
    if gp.selection.tournament.lex_pressure
        lex_inf = 'True';
    else
        lex_inf = 'False';
    end
    
    disp(' ');
    disp('-------------------------------------------------------------------------');
    disp('GPTIPS 2F');
    disp('Symbolic data mining platform for MATLAB evolved');
    disp('Copyright (C) Dominic Searson* 2009-2015, Aleksei Tepljakov& 2017-');
    disp(' ');
    disp('Contact: * - searson@gmail.com, & - alex@starspirals.net');
    disp(' ');
    
    disp('This program is free software: you can redistribute it and/or modify');
    disp('it under the terms of the GNU General Public License as published by');
    disp('the Free Software Foundation, either version 3 of the License, or');
    disp('(at your option) any later version.');
    disp(' ');
    disp('This program is distributed in the hope that it will be useful,');
    disp('but WITHOUT ANY WARRANTY; without even the implied warranty of');
    disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
    disp('GNU General Public License for more details: http://www.gnu.org/licenses');
    disp(' ');
    disp('-------------------------------------------------------------------------- ');
    disp(' ');
    disp('Run parameters');
    disp('--------------');
    disp(['Population size:         ' int2str(gp.runcontrol.pop_size)]);
    disp(['Number of generations:   ' int2str(gp.runcontrol.num_gen)]);
    disp(['Number of runs:          ' int2str(gp.runcontrol.runs)]);
    
    if gp.runcontrol.parallel.auto
        disp('Parallel mode :          auto ');
    elseif gp.runcontrol.parallel.enable
        disp('Parallel mode :          manual ');
    else
        disp('Parallel mode :          off ');
    end
    
    if gp.selection.tournament.p_pareto == 0
        disp('Tournament type:         regular');
    else
        disp(['Tournament type:         Pareto (probability = ' ...
            num2str(gp.selection.tournament.p_pareto) ')']);
    end
    
    disp(['Tournament size:         ' int2str(gp.selection.tournament.size)]);
    disp(['Elite fraction:          ' num2str(gp.selection.elite_fraction)]);
    
    if gp.runcontrol.usecache
        disp('Fitness cache:           enabled');
    else
        disp('Fitness cache:           disabled');
    end
    
    disp(['Lexicographic selection: ' lex_inf]);
    disp(['Max tree depth:          ' int2str(gp.treedef.max_depth)]);
    disp(['Max nodes per tree:      ' int2str(gp.treedef.max_nodes)]);
    disp(['Using function set:     ' fns]);
    
    % If ADFs are used, display the list here
    % TODO: maybe overwhelming? in that case, limit to, e.g., 10-15 entries
    if gp.nodes.adf.use
        disp('List of ADF functions:');
        for k=1:length(gp.nodes.adf.name)
           disp(['   ' gp.nodes.adf.name{k} ': ' gp.nodes.adf.eval{k}]); 
        end
    end
    
    % If evolutionary rules are applied, list them at this point
    if gp.evolution.rules.use
        disp('Using evolutionary rules:');
        for k=1:length(gp.evolution.rules.sets)
            the_rule = gp.evolution.rules.sets{k};
            disp(['   ' func2str(the_rule{1})]);
        end
    end
    
    disp(['Number of inputs:        ' int2str(gp.nodes.inputs.num_inp)]);
    
    if gp.genes.multigene
        disp(['Max genes:               ' int2str(gp.genes.max_genes)]);
    end
    
    if ~gp.nodes.const.p_ERC
        disp('Using no constants');
    else
        disp(['Constants range:         [' num2str(gp.nodes.const.range) ']']);
    end
    
    if gp.nodes.pconst.p_PRC < eps
        disp('Using no preset constants');
    else
        disp(['Preset Constants array:  [' sprintf('%.2f ', gp.nodes.pconst.set) ']']);
    end
    
    if gp.fitness.complexityMeasure
        disp('Complexity measure:      expressional');
    else
        disp('Complexity measure:      node count');
    end
    
    disp(['Fitness function:        ' func2str(gp.fitness.fitfun) '.m']);
    disp(' ');
end

%parallel computing initialisation
gp = gpinitparallel(gp);

%log run start time
gp.info.startTime = datestr(now,0);

%set elapsed run time count to zero seconds
gp.state.runTimeElapsed = 0;

function gp=procadfs(gp)
%PROCADFS Processes ADF entries (if any), assigns names to these functions
%and adds them to the list of function candidates along with corresponding
%"active" parameters, if present

% Perform additional check: if there are no ADFs, then there is nothing to do.
exprs = gp.nodes.adf.expr;
numExprs = length(exprs);
if isempty(exprs)
    gp.nodes.adf.use = false;
    return;
end
gp.nodes.adf.use = true;

% Unless names are specified, generate the names automatically
adfname = 'adf';
use_auto_adf_names = false;
 
if ~isfield(gp.nodes.adf, 'name') || isempty(gp.nodes.adf.name)
    gp.nodes.adf.name = cell(1,numExprs);
    use_auto_adf_names = true;
end  

% Process the expressions
gp.nodes.adf.eval = cell(1,numExprs);
gp.nodes.adf.seed = cell(1,numExprs);
for k=1:numExprs
    [fcn, ~, ~, seed] = parseadf(exprs{k});
    if use_auto_adf_names, gp.nodes.adf.name{k} = [adfname num2str(k)]; end
    gp.nodes.adf.eval{k} = fcn;
    gp.nodes.adf.seed{k} = seed;
end

% We now add the ADFs as regular functions to be used in GP

% Check for active ADFs
if ~isfield(gp.nodes.adf,'active') || isempty(gp.nodes.adf.active)
    gp.nodes.adf.active = ones(1,length(gp.nodes.adf.name));
end

% Add ADFs to the current function list
gp.nodes.functions.name = ...
    [gp.nodes.functions.name(:)' gp.nodes.adf.name(:)'];

% Enumerate all expressions found in ADFs, because their symbolic names
% must be known to MATLAB when converting to symbolic expressions
funcs = {};
for k=1:numExprs
    out = regexp(exprs{k}, '([A-Za-z0-9\_]+)( *)?\(', 'tokens');
    for l=1:numel(out)
        part = out{l};
        funcs{end+1} = part{1}; %#ok - no need to parse twice
    end     
end

% Finally, store this data for later use in gpsym
gp.nodes.adf.syms = unique(funcs);


function gp=procfuncnodes(gp)
%PROCFUNCNODES Process required function node information prior to a run.

% Account for generated ADF functions' arity. We need the corresponding
% handles defined in the current workspace for that reason
if gp.nodes.adf.use, assignadf(gp); end

%loop through function nodes and generate arity list
for i = 1:length(gp.nodes.functions.name)
    
    % nargin does not detect anonymous functions when string is supplied.
    % We'll have to check for these functions manually.
    if exist(gp.nodes.functions.name{i},'var') % Variable in the workspace
        arity = nargin(eval(gp.nodes.functions.name{i}));
    else
        arity = nargin(gp.nodes.functions.name{i});
    end
    
    %some functions have a variable number of input arguments (e.g. rand)
    %In this case generate an error message and exit
    if arity == -1
        error(['The function ' gp.nodes.functions.name{i} ...
            ' may not be used (directly) as a function node because it has a variable number of arguments.']);
    end
    
    gp.nodes.functions.arity(i) = arity;
end

if ~isfield(gp.nodes.functions,'active') || isempty(gp.nodes.functions.active)
    gp.nodes.functions.active = ones(1,length(gp.nodes.functions.name)-length(gp.nodes.adf.name));
end

% If ADFs are used, add their active functions
if gp.nodes.adf.use
    gp.nodes.functions.active = [gp.nodes.functions.active gp.nodes.adf.active];
end

gp.nodes.functions.active = logical(gp.nodes.functions.active);

%check max number of allowed functions not exceeded
gp.nodes.functions.num_active = numel(find(gp.nodes.functions.active));
if gp.nodes.functions.num_active > 22
    error('Maximum number of active functions allowed is 22');
end

% Get number of active ADFs, if present
if gp.nodes.adf.use
   gp.nodes.adf.num_active = numel(find(gp.nodes.adf.active)); 
end

%Generate single char Active Function IDentifiers (afid)(a->z excluding
%x,e,i,j) to stand in for function names whilst processing expressions.
%Exclusions are because 'x' is reserved for input nodes, 'e' is used for
%expressing numbers in standard form by Matlab and, by default, 'i' and 'j'
%represent sqrt(-1).
charnum = 96; skip = 0;
for i=1:gp.nodes.functions.num_active
    while true      %e                          %i                    %j                         %x
        if (charnum+i+skip)==101 || (charnum+i+skip)==105 || (charnum+i+skip)==106 || (charnum+i+skip)==120
            skip = skip + 1;
        else
            break
        end
        
    end
    afid(i) = char(charnum+i+skip);
end

%extract upper case active function names for later use
gp.nodes.functions.afid = afid;

if numel(gp.nodes.functions.name) ~= numel(gp.nodes.functions.active)
    error('There must be the same number of entries in gp.nodes.functions.name and gp.nodes.functions.active. Check your config file.');
end

% We need to populate the seeds used in seeding, if ADFs are defined
% ADFs always appear at the end, so we can get the corresponding characters
% easily from there. However, we also need to remove unnecessary seeds from
% the initial seed pool.
if gp.nodes.adf.use
    use_seeds = {};
    the_seeds = cell(gp.nodes.adf.num_active,1);
    [the_seeds{:}] = deal(gp.nodes.adf.seed{logical(gp.nodes.adf.active)});
    seedind = gp.nodes.functions.num_active - gp.nodes.adf.num_active;
    seedstr = '';
    for k=1:length(the_seeds)
        use_seeds{k} = [gp.nodes.functions.afid(seedind+k) the_seeds{k}];
        seedstr = [seedstr gp.nodes.functions.afid(seedind+k)];
    end
    gp.nodes.adf.use_seeds = use_seeds;
    gp.nodes.adf.seed_str = seedstr;
end

temp = cell(gp.nodes.functions.num_active,1);
[temp{:}] = deal(gp.nodes.functions.name{gp.nodes.functions.active});
[gp.nodes.functions.active_name_UC] = upper(temp);

%generate index locators for arity >0 and arity == 0 active functions. The
%treegen function needs this info later for identifying which functions are
%terminal and which are internal nodes.
active_ar = (gp.nodes.functions.arity(gp.nodes.functions.active));
fun_argt0 = active_ar > 0;
fun_areq0 =~ fun_argt0;

gp.nodes.functions.afid_argt0 = gp.nodes.functions.afid(fun_argt0); %functions with arity > 0
gp.nodes.functions.afid_areq0 = gp.nodes.functions.afid(fun_areq0); %functions with arity == 0
gp.nodes.functions.arity_argt0 = active_ar(fun_argt0);

gp.nodes.functions.fun_lengthargt0 = numel(gp.nodes.functions.afid_argt0);
gp.nodes.functions.fun_lengthareq0 = numel(gp.nodes.functions.afid_areq0);

function gp = initcache(gp)
%INITCACHE Sets up fitness cache.
gp.fitness.cache = containers.Map('keytype','uint32','valuetype','any');
