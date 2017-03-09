% GPTIPS2
% Symbolic Data Mining for MATLAB evolved
% (c) Aleksei Tepljakov 2017-...
% (c) Dominic Searson 2009-2015
%
% Files
%   bootsample                    - Get an index vector to sample with replacement a data matrix X.
%   comparemodelsREC              - Graphical REC performance curves for between 1 and 5 multigene models.
%   crossover                     - Sub-tree crossover of encoded tree expressions to produce 2 new ones.
%   displaystats                  - Displays run stats periodically.
%   drawtrees                     - Draws the tree structure(s) of an individual in a web browser.
%   evalfitness                   - Calls the user specified fitness function.
%   evalfitness_par               - Calls the user specified fitness function (parallel version).
%   extract                       - Extract a subtree from an encoded tree expression.
%   genebrowser                   - Visually analyse unique genes in a population and identify horizontal bloat.
%   genefilter                    - Removes highly correlated genes from a unique GENES struct.
%   genes2gpmodel                 - Create a data structure representing a multigene symbolic regression model from the specified gene list.
%   getcomplexity                 - Returns the expressional complexity of an encoded tree or a cell array of trees.
%   getdepth                      - Returns the tree depth of an encoded tree expression.
%   getnumnodes                   - Returns the number of nodes in an encoded tree expression or the total node count for a cell array of expressions.
%   gp_2d_mesh                    - Creates new training matrix containing pairwise values of all x1 and x2.
%   gpcheck                       - Perform pre-run error checks.
%   gpdefaults                    - Initialises the GPTIPS struct by creating default parameter values.
%   gpfinalise                    - Finalises a run.
%   gpinit                        - Initialises a run.
%   gpinitparallel                - Initialise the Parallel Computing Toolbox.
%   gpmodelfilter                 - Object to filter a population of multigene symbolic regression models.
%   gpmodelreport                 - Generate an HTML report on the specified multigene regression model.
%   gpmodelvars                   - Display the frequency of input variables present in the specified model.
%   gppopvars                     - Display frequency of the input variables present in models in the population.
%   gppretty                      - Simplify and prettify a multigene symbolic regression model.
%   gprandom                      - Sets random number generator seed according to system clock or user seed.
%   gpreformat                    - Reformats encoded trees so that the Symbolic Math toolbox can process them properly.
%   gpsimplify                    - Simplify SYM expressions in a less glitchy way than SIMPLIFY or SIMPLE.
%   gpterminate                   - Check for early termination of run.
%   gptic                         - Updates the running time of this run.
%   gptoc                         - Updates the running time of this run.
%   gptoolboxcheck                - Checks if certain toolboxes are installed and licensed.
%   gptreestructure               - Create cell array containing tree structure connectivity and label information for an encoded tree expression.
%   initbuild                     - Generate an initial population of GP individuals.
%   kogene                        - Knock out genes from a cell array of tree expressions.
%   mergegp                       - Merges two GP population structs into a new one.
%   mutate                        - Mutate an encoded symbolic tree expression.
%   ndfsort_rank1                 - Fast non dominated sorting algorithm for 2 objectives only - returns only rank 1 solutions. 
%   paretoreport                  - Generate an HTML performance/complexity report on the Pareto front of the population.
%   picknode                      - Select a node (or nodes) of specified type from an encoded GP expression and return its position. 
%   popbrowser                    - Visually browse complexity and performance characteristics of a population.
%   popbuild                      - Build next population of individuals.
%   pref2inf                      - Recursively extract arguments from a prefix expression and convert to infix where possible.
%   regressionErrorCharacteristic - Generates REC curve data using actual and predicted output vectors.
%   rungp                         - Runs GPTIPS 2 using the specified configuration file.
%   runtree                       - Run the fitness function on an individual in the current population.
%   scangenes                     - Scan a single multigene individual for all input variables and return a frequency vector.
%   selection                     - Selects an individual from the current population.
%   standaloneModelStats          - Compute model performance stats for actual and predicted values.
%   summary                       - Plots basic summary information from a run.
%   tree2evalstr                  - Converts encoded tree expressions into math expressions that MATLAB can evaluate directly.
%   treegen                       - Generate a new encoded GP tree expression.
%   uniquegenes                   - Returns a GENES structure containing the unique genes in a population.
%   updatestats                   - Update run statistics.
