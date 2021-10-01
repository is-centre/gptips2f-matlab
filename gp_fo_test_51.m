function gp = gp_fo_test_51(gp)

%% The function under test

%g = fotf('s^1.2+3s^0.4+5', 's^1.6+10s^1.2+35s^0.8+50s^0.4+24')

%g = fotf('s^3+3s+5', 's^4+10s^3+35s^2+50s+24'); g = oustapp(g) % This
%system is not working well with this approach!

% New IO system to test ...
%g = fotf('11s^{3}+179s^{2}+757s+717', 's^{4}+18s^{3}+104s^{2}+222s+135'); g=oustapp(g);

% Repeated poles
g = fotf('1', 's^1.2+5s^0.9+9s^0.6+7s^0.3+2')

%% Generate the identification data
% In addition, define the global dt_conv parameter
global dt_conv;
dt = 0.01;
dt_conv = dt;

% Ref signals: C change at the end required by nonlinearities in U,
% otherwise it could potentially remain undetected by the identification
% algorithm
st = refgen('prbs', 1, 50, 'sin', 1.5, 20, 'c', 1, 25, 'c', 1.5, 25, 'c', 0.5, 25, 'Ts', dt);
stt = refgen('c', 1, 15, 'sin', 2.5, 10, 'prbs', [1 1 0 15], 50,  'Ts', dt);

%plot(st);

u = st.u;
u1 = stt.u;
nonlin = true;

% Try Input nonlinearities
if nonlin
    u = st.u.^2;
    u1 = stt.u.^2;
end

% Training
y = lsim(g, u, st.t);

% Test
y1 = lsim(g,u1, stt.t);

%noiseVal = 1e-2; % For IO
noiseVal = 1e-5; % For FO

% Add some noise (or not)
y = y + (noiseVal*(-(2*randn(size(st.u)))));
y1 = y1 + (noiseVal*(-(2*randn(size(stt.u)))));

gp.userdata.xtrain = [st.t st.u];
gp.userdata.ytrain = y;

% Validation
gp.userdata.xval = [stt.t stt.u];
gp.userdata.yval = y1;
gp.userdata.user_fcn = @regressmulti_fitfun_validate;

% Test
gp.userdata.xtest = [stt.t stt.u];
gp.userdata.ytest = y1;

gp.userdata.name = 'Dynamic system identification test with ADFs';

%% Functions
gp.nodes.functions.name = {'times','minus','plus','square','sqrt'};

% Assign names to ADFs
gp.nodes.adf.name = {'simp_f', 'simp_f_c1', 'simp_f_c2', 'simp_exp'};

% ADFs for repeated poles
gp.nodes.adf.expr = {'convimp($1, fian(times(?1,$2),fila(#1)))', ...
    'convimp($1, convself(fian(times(?1,$2),fila(#1))))', ...
    'convimp($1, convself2(fian(times(?1,$2),fila(#1))))', ...
   'convimp($1, exp(times(?1,$2)))'};

gp.nodes.pconst.set = [0.1 0.2 0.3 0.4 0.5];
gp.nodes.pconst.p_PRC = 0.25;
gp.nodes.adf.p_gen = 0.25;

% Enforce structure on mutation
gp.nodes.adf.arg_force = true;

%run control parameters
gp.runcontrol.pop_size = 500;                     
gp.runcontrol.num_gen = 100;				                                                 

% Initial setting
gp.runcontrol.timeout = 360;

% Selection
gp.selection.tournament.size = 50;
gp.selection.tournament.p_pareto = 0.3;
gp.selection.elite_fraction = 0.2;      

% Termination criterion
gp.fitness.terminate = true;
gp.fitness.terminate_value = 1e-4;

% Trees
gp.treedef.max_depth = 3;
gp.treedef.max_mutate_depth = 3;
gp.treedef.max_nodes = 12;

% Mutation
%gp.operators.mutation.mutate_par = [0.8 0.05 0.1 0 0.05 0];

gp.nodes.inputs.names = {'t', 'u'};

%genes
gp.genes.max_genes = 5;

gp.operators.mutation.p_mutate = 0.24;    
gp.operators.crossover.p_cross = 0.71;    
gp.operators.directrepro.p_direct = 0.05;

% Constant nodes
gp.nodes.const.range = [-10 0];
gp.nodes.const.p_ERC = 0.4;
gp.nodes.const.p_int = 0.3;

