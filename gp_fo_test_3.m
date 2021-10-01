function gp = gp_fo_test_3(gp)
%Try symbolic regression for FO dynamic system identification

%run control parameters
gp.runcontrol.pop_size = 100;                     
gp.runcontrol.num_gen = 100;				                                                 
			
%selection
gp.selection.tournament.size = 6;

%termination
gp.fitness.terminate = true;
gp.fitness.terminate_value = 0.003;

% The function under test
g = fotf('s^1.2+3s^0.4+5', 's^1.6+10s^1.2+35s^0.8+50s^0.4+24')
%g = fotf('s^2+3s+5', 's^4+10s^3+35s^2+50s+24')

% Generate the identification data
dt = 0.01;
st = refgen('prbs', 1, 50, 'Ts', dt);
u = st.u;
t = st.t;
y = lsim(g, u, t);

% Add some noise
y = y + (1e-6*(-(2*randn(size(u)))));

% In addition, define the global dt_conv parameter
global dt_conv;
dt_conv = dt;

gp.userdata.xtrain = [t u];
gp.userdata.ytrain = y;
gp.userdata.xtest = [t u];
gp.userdata.ytest = y;
gp.userdata.name = 'Dynamic system identification test';

%genes
gp.genes.max_genes = 5;                   

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
    'sin','cos','exp','add3','mult3', ...
    'convimp','convself', 'convself2',...
    'f_alpha_0_1','f_alpha_0_4'};