%% This script is used to generate Pade-like approximations of the MLF
%  function of the E_alpha,alpha(-t) kind for several fixed values of
%  alpha. The application of these approximations is as follows: to
%  allow to create candidate functions for time domain identification
%  that are efficiently computed. Although the precision is lost, at least
%  initial evaluation of the candidate models will be sufficiently fast.

clear
clc

%% We use the MLF function implementation by Podlubny as basis

% Here we set the parameters of the time-domain region of interest

% Magnitudes & number of points
t_mag_l = 10^-5;
t_mag_h = 10^5;
mlf_prec = 10;

% Time step
dt = 0.01;

% Number of points
NumPts = floor((t_mag_h - t_mag_l)/dt);

%% The parameters of approximants
m = 1:3;
n = 1:7;
alpha = 0.05:0.05:0.95;

%% Generate the time vector
t = linspace(t_mag_l, t_mag_h, NumPts);

%% Gentlemen, start your engines!
testNum = 1;
nowT = tic;
for italpha = alpha
   
    disp(['Starting tests for alpha=' num2str(italpha) '...']);
    y = mlf(italpha,italpha,-t,mlf_prec);
    
    for itm = m
        for itn = n
            
            rzer = 1:itm
            rpol = (1+itm):(itm+itn)
            p = ones(itm+itn,1);
            
            % Nonlinear least-squares using Levenberg-Marquardt algorithm
            f1 = @(z) y - polyval(z(rzer),-t)./polyval(z(rpol),-t);
            
            % Approximated function
            f_ = @(z,c) polyval(z(rzer),-c) ./ polyval(z(rpol), -c);
            
            xo = lsqnonlin(f1, p, [], [], ...
                optimset('Display', 'iter', 'MaxFunEvals', 20000, ...
                         'MaxIter', 10000, ...
                         'Algorithm', 'levenberg-marquardt'));
                     
            % Do all sorts of error computation here
            err = y-f_(xo,t);               % Residuals as is
            err_abs_sum = sum(abs(err));    % Sum of absolute error
            err_rel = abs(err ./ y)*100;    % Relative error [%]
            err_rel_max = max(err_rel);     % Maximum relative error
            
            % Store EVERYTHING
            the_results = struct;
            
            the_results.t = t;
            the_results.y = y;
            the_results.alpha = italpha;
            the_results.b = xo(rzer);
            the_results.a = xo(rpol);
            the_results.err = err;
            the_results.err_abs_sum = err_abs_sum;
            the_results.err_rel = err_rel;
            the_results.err_rel_max = err_rel_max;
            
            results{testNum} = the_results;
            testNum = testNum + 1;
        end
    end
    
    disp('Done.');
end
disp('The tests are over. Total running time [in seconds] below');
toc(nowT)

disp('Saving the file as well.');
save mlf_approx_results.mat results