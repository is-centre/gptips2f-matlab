%% INITIALIZATION

% Function names and file name patterns
fun_names = {'int_01_fun', 'frac_01_fun', 'frac_02_fun'};
file_patt = '_test_%04d'; % Added to fun name to save the mat file
table_patt = '_table'; % The variable name for corresponding table
folder_path = 'exp_results_v1';

% Random number generator initialization to ensure repeatability
rng(2332);

% Define number of tests; gp config is saved for each test
numTests = 100;

%% TESTS

% Create the folder
mkdir(folder_path);

% Conduct the tests for all the functions
for k=1:length(fun_names)
    % Function handle
    hand = str2func(fun_names{k});
    for l=1:numTests
        gp = rungp(hand);
        fileName = [fun_names{k} sprintf(file_patt, l)];
        filePath = [folder_path '\' fileName '.mat'];
        save(filePath, 'gp');
    end
end

%% STATISTICS

% Load all files one by one and get statistical data into a single
% table for each function

% Define the global dt_conv: in case this script is used later to only
% process the available data
global dt_conv;
dt_conv = 0.01;

for k=1:length(fun_names)
    
    T = table;
    
    % Count the times 'correct' nodes are selected among ALL nodes
    % NB!!! This is specific to the test, so be careful if you want to
    % reuse this piece of code later, do so with caution
    all_nodes = {'frac', 'cint'};
    switch k
        case 1
            corr_node = 2; % Reference to item (index) in cell array
        case 2
            corr_node = 1;
        case 3
            corr_node = 1;
        otherwise
            error('Check your code!');
    end
    
    
    % We use VALBEST parameter since validation allows to combat
    % overfitting and is therefore more important than training fitness
    for l=1:numTests
        fileName = [fun_names{k} sprintf(file_patt, l)];
        filePath = [folder_path '\' fileName '.mat'];
        load(filePath)
        
        % Get all the interesting parameters out of this test.
        % We'll use valbest because it describes
        gp_valfit = gp.results.valbest.valfitness;
        gp_comp = gp.results.valbest.complexity;
        gp_nodecnt = gp.results.valbest.nodecount;
        
        % Get residuals information
        [er,st] = gpresids(gp,'valbest');
        gp_resid_max_ac = max(st.RErr/st.ConfBnd);
 
        % Go through eval_individual strings and find the nodes
        
        % Row vector of counts where the indices are the same as in
        % all_nodes cell array
        the_node_counts = zeros(1,length(all_nodes));
        for m=1:length(gp.results.valbest.eval_individual)
            for n=1:length(all_nodes)
                this_cnt = strfind(gp.results.valbest.eval_individual{m}, ...
                    all_nodes{n});
                the_node_counts(n) = the_node_counts(n) + length(this_cnt);
            end
        end
        
        % Find the statistics
        gp_corr_nodes = the_node_counts(corr_node);
        gp_corr_node_perc = gp_corr_nodes / sum(the_node_counts) * 100;
        
        % Add the corresponding row to table
        varNames = {'TestNo','Fitness','MaxResidAC','SimilarNodes','SimilarNodesPerc','Complexity','NodeCount'};
        varVals = {l, gp_valfit, gp_resid_max_ac, gp_corr_nodes, gp_corr_node_perc, gp_comp, gp_nodecnt};
        T = [T; cell2table(varVals, 'VariableNames', varNames)];
        
    end
    
    % Sort rows by fitness
    % T = sortrows(T,'Fitness','ascend');
    
    % Alternatively, by MaxResidAC's
    T = sortrows(T,'MaxResidAC','ascend');
    
    fileName = [fun_names{k} table_patt];
    filePath = [folder_path '\' fileName '.mat'];
    save(filePath, 'T');
    
end