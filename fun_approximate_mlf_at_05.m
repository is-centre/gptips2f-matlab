function gp = fun_approximate_mlf_at_05(gp)
%Try symbolic regression for FO dynamic system identification

%run control parameters
gp.runcontrol.pop_size = 100;                     
gp.runcontrol.num_gen = 100;				                                                 
			
%selection
gp.selection.tournament.size = 6;

%termination
gp.fitness.terminate = true;
gp.fitness.terminate_value = 1e-12;

% The data!
disp('Start generating points for MLF function ... ');
alpha = 0.5;
x = linspace(10^(-5),10^5,10^6);
prec = 10;

% Get the response from the original function
y = mlf(alpha,alpha,-x,prec);

disp('Done.');
disp(' ');

x = colv(x);
y = colv(y);

gp.userdata.xtrain = x;
gp.userdata.ytrain = y;
gp.userdata.xtest = x;
gp.userdata.ytest = y;
gp.userdata.name = ['MLF identification for alpha=' num2str(alpha)];

%genes
gp.genes.max_genes = 3;                   

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','sqrt','square', ...
    'sin','cos','exp','add3','mult3'};