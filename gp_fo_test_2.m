function gp = gp_fo_test_2(gp)
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
loadsets

% Generate the identification data
dt = 0.01;
st = refgen('prbs', 1, 50, 'Ts', dt);
u = st.u;
t = st.t;
y = lsim(G2, u, t);

% In addition, define the global dt_conv parameter
global dt_conv;
dt_conv = dt;

gp.userdata.xtrain = [t u];
gp.userdata.ytrain = y;
gp.userdata.xtest = [t u];
gp.userdata.ytest = y;
gp.userdata.name = 'Dynamic system identification test';

%genes
gp.genes.max_genes = 3;                   

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
    'sin','cos','exp','add3','mult3', ...
    'convimp','f_alpha_0_1','f_alpha_0_25'};