% Newly implemented ADF test. In fact, we are going to do many tests with
% different ideas and such ... let's see how well we can identify our
% systems!

%% The identification part
disp('Running a more difficult FO TF identification task');
disp('using Symbolic Regression through GPTIPS2 package.');
disp(' ');
gp=rungp(@gp_fo_test_52);

%% Summary
disp('Next, plot summary information of run using:');
disp('>>summary(gp)');
disp(' ');
summary(gp,false);

%% runtree
disp('Run the best model of the run on the fitness function using:');
disp('>>runtree(gp,''best'');');
disp(' ');
runtree(gp,'best');

%% Gppretty
disp('Best simplified model');
gppretty(gp,'best');

disp(['This model has a tree depth of ' int2str( getdepth(gp.results.best.individual{1}))]);
disp(['It was found at generation ' int2str(gp.results.best.foundatgen)]);
disp(['and has fitness ' num2str(gp.results.best.fitness)]);

%% Resids
disp('Showing residuals for training and validation data');
gpresids(gp,'best');
gpresids(gp,'valbest');

%% Tree
disp(' ');
disp('Next, visualise the tree structure of the best model of the run:');
disp('>>drawtrees(gp,''best'');');
disp(' ');
drawtrees(gp,'best');