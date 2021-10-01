function gp = gp_fo_test_4(gp)
%Try symbolic regression for FO dynamic system identification

%run control parameters
gp.runcontrol.pop_size = 250;                     
gp.runcontrol.num_gen = 100;				                                                 

% Two minutes is all that's needed (for our initial testing)
gp.runcontrol.timeout = 120;

%selection
gp.selection.tournament.size = 15;
gp.selection.tournament.p_pareto = 0.3; % These are new
gp.selection.elite_fraction = 0.3;      % parameters

%termination
gp.fitness.terminate = true;
gp.fitness.terminate_value = 1e-4;

% Trees
gp.treedef.max_depth = 4;
gp.treedef.max_mutate_depth = 4;
gp.treedef.max_nodes = 12;

% Mutation
% gp.operators.mutation.mutate_par = [0.8 0.05 0.1 0 0.05 0];

gp.nodes.inputs.names = {'u', 't'};

%genes
gp.genes.max_genes = 5;

gp.operators.mutation.p_mutate = 0.24;    
gp.operators.crossover.p_cross = 0.74;    
gp.operators.directrepro.p_direct = 0.02;

% Constant nodes
gp.nodes.const.range = [-10 10];
gp.nodes.const.p_ERC = 0.4;
gp.nodes.const.p_int = 0.3;

%define function nodes
% gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
%     'exp','add3','mult3','convimp','convself','convself2','fian_01', ...
%     'fian_04','fian_05', 'sin', 'cos'};

%  gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
%      'exp','add3','mult3','convimp','convself','convself2','fian_01', ...
%      'fian_03','fian_04'};
 
  %gp.nodes.functions.name = {'times','minus','plus', 'square', 'sqrt', ...
  %   'exp','convimp','convself','convself2','fian_05'};
 
 % Yet another version of this with two equivalent nodes to describe
 % dynamics
 gp.nodes.functions.name = {'times','minus','plus', 'square', 'sqrt', ...
     'convimp','convself','convself2','co_ut','cof_ut_05'};

  %gp.nodes.functions.name = {'times','minus','plus', 'exp','convimp', 'fian_05'};
  
  % Essentially, this is cheating. Let machine
  % learning figure stuff out on its own.
  %gp.nodes.functions.name = {'times','minus','plus','exp','sqrt','square', ...
  %    'cof_ut_01', 'cof_ut_03', 'cof_ut_04', 'cof_ut_05', 'co_ut', ...
  %    'convself', 'convself2'};

%gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
%    'exp','add3','mult3','convimp'};

% The function under test
%g = fotf('s^1.2+3s^0.4+5', 's^1.6+10s^1.2+35s^0.8+50s^0.4+24')

%g = fotf('s^3+3s+5', 's^4+10s^3+35s^2+50s+24'); g = oustapp(g) % This
%system is not working well with this approach!

% New IO system to test ...
g = fotf('11s^{3}+179s^{2}+757s+717', 's^{4}+18s^{3}+104s^{2}+222s+135'); g=oustapp(g);

% Repeated poles
%g = fotf('1', 's^1.2+5s^0.9+9s^0.6+7s^0.3+2')

% Generate the identification data
dt = 0.01;
%st = refgen('prbs', 1, 50, 'Ts', dt);
st = refgen('prbs', 1, 50, 'sin', 1.5, 20, 'c', 0.5, 25, 'Ts', dt);
%st = refgen('prbs', 5, 150, 'Ts', dt); 
stt = refgen('sin', 2.5, 10, 'prbs', [1 1 0 15], 50, 'c', 0.75, 15, 'Ts', dt);
%stv = refgen('prbs', [1 1 0 27], 100, 'Ts', dt);
u = st.u;
t = st.t;
y = lsim(g, u, t);

% figure; subplot(211); plot(t,y); subplot(212); plot(t,st.u); drawnow;


% Training
y1 = lsim(g,stt.u, stt.t);
%y2 = lsim(g,stv.u, stv.t);

% Add some noise (or not)
%y = y + (1e-3*(-(2*randn(size(u)))));
%y1 = y1 + (1e-3*(-(2*randn(size(stt.u)))));
%y2 = y2 + (1e-5*(-(2*randn(size(stt.u)))));

%figure; plot(t,y); drawnow;

% In addition, define the global dt_conv parameter
global dt_conv;
dt_conv = dt;

gp.userdata.xtrain = [t u];
gp.userdata.ytrain = y;

% Validation
%gp.userdata.xval = [stv.t stv.u];
%gp.userdata.yval = y2;

% Test
gp.userdata.xtest = [stt.t stt.u];
gp.userdata.ytest = y1;

gp.userdata.name = 'Dynamic system identification test';