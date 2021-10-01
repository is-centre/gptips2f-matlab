function gp = gp_fo_test_1(gp)
%Try symbolic regression for FO dynamic system identification

%run control parameters
gp.runcontrol.pop_size = 100;                     
gp.runcontrol.num_gen = 100;				                                                 
			
%selection
gp.selection.tournament.size = 6;

%termination
gp.fitness.terminate = true;
gp.fitness.terminate_value = 0.003;

%load in the raw x and y data
load gpid1_test

% In addition, define the global dt_conv parameter
global dt_conv;
dt_conv = gpid1.dt;

gp.userdata.xtrain = [gpid1.t gpid1.u];
gp.userdata.ytrain = gpid1.y;
gp.userdata.xtest = [gpid1.t gpid1.u];
gp.userdata.ytest = gpid1.y;
gp.userdata.name = 'Dynamic system identification test';

%genes
gp.genes.max_genes = 3;                   

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
    'sin','cos','exp','add3','mult3', ...
    'convimp','f_alpha_0_1','f_alpha_0_25'};