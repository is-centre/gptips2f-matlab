clc;
disp('Let''s try actually approximating the MLF(a,a,x) for a=0.5');
disp('using Symbolic Regression through GPTIPS2 package.');
disp('Press a key to continue');
disp(' ');
pause;
gp=rungp(@fun_approximate_mlf_at_05);

disp('Next, plot summary information of run using:');
disp('>>summary(gp)');
disp('Press a key to continue');
disp(' ');
pause;
summary(gp,false);

disp('Run the best model of the run on the fitness function using:');
disp('>>runtree(gp,''best'');');
disp('Press a key to continue');
disp(' ');
pause;
runtree(gp,'best');

disp(' ');
disp('The best model of the run is stored in the field:');
disp('gp.results.best.eval_individual{1} :');
disp(' ');
disp( gp.results.best.eval_individual{1});
disp(' ');

disp(['This model has a tree depth of ' int2str( getdepth(gp.results.best.individual{1}))]);
disp(['It was found at generation ' int2str(gp.results.best.foundatgen)]);
disp(['and has fitness ' num2str(gp.results.best.fitness)]);

%If Symbolic Math toolbox is present
if gp.info.toolbox.symbolic
    
    disp(' ');
    disp('Using the symbolic math toolbox simplified versions of this');
    disp('expression can be found: ')
    disp('E.g. using the the GPPRETTY command on the best model: ');
    disp('>>gppretty(gp,''best'') ');
    disp('Press a key to continue');
    disp(' ');
    pause;
    gppretty(gp,'best');
    disp(' ');
    
end

disp(' ');
disp('Next, visualise the tree structure of the best model of the run:');
disp('>>drawtrees(gp,''best'');');
disp('Press a key to continue');
disp(' ');
pause;
drawtrees(gp,'best');