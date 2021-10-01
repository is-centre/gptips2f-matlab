% Here we simply test the MLF function for solving some FODE's

% We're going to use the following equation
%             1/a 
% G(s) = -------------
%        s^alpha + b/a

format long

% Define parameters
a = 1;
b = 0.5;
alpha = 0.1;

% Time
h = 0.01;
t_max = 100;
t = 0:h:t_max;

%% Test the MLF approximations!

% Less time: up to 10s
t1 = t;
t1(t1>10)=[];

tm = tic; p1 = mlf(alpha, alpha, -t1, 10); toc(tm)
tm = tic; p2 = mlf_gpa_a(alpha, t1); toc(tm)
%p3 = ml(-t,alpha);

figure;
subplot(211);
plot(t1, p1); hold on; plot(t1, p2); %plot(p3);
legend('mlf', 'gpa'); %, 'grap');
grid;
title(['ML function approximations, alpha: ' num2str(alpha)]);
subplot(212);
aerr = sum(abs(p1-p2));
strrep(num2str(aerr, 10), '.', ',')
plot(t1, p1-p2);
grid;
title(['Residuals; sum of absolute error: ' num2str(aerr)]);

%% Define and analyse fotf
G = 1/a / (fotf('s')^alpha + b/a);
u = ones(size(t));
u(t>(t_max/2)) = 0.5;
y = lsim(G,u,t);
y_y = impulse(G,t);

%% Will Y be computed the same here?
y = (t(2)-t(1)) .* conv(u, y_y);
y = y(1:numel(t));

%% Now, do the same using convolution/MLF function
tm = tic;
y_imp = fo_impulse_f1(t,b,alpha);
y_mlf = (t(2)-t(1)) .* conv(u, y_imp);
y_mlf = y_mlf(1:numel(t));
tm = toc(tm);
disp(['MLF computation: ' num2str(tm)]);

%% And its approximation
tm = tic;
y_imp1 = fo_impulse_f1_gpa(t,b,alpha);
y_mlf1 = (t(2)-t(1)) .* conv(u, y_imp1);
y_mlf1 = y_mlf1(1:numel(t));
tm = toc(tm);
disp(['GPA-ML computation: ' num2str(tm)]);

%% And one more
tm = tic;
y_imp2 = fo_impulse_f1_ml(t,b,alpha);
y_mlf2 = (t(2)-t(1)) .* conv(u, y_imp2);
y_mlf2 = y_mlf2(1:numel(t));
tm = toc(tm);
disp(['ML computation: ' num2str(tm)]);

%% Compare ALL
figure;
plot(y_mlf); hold on;
plot(y_mlf1);
plot(y_mlf2);
plot(y);
legend('mlf', 'gpa', 'gra', 'orig');

%% Compare fixed MLF1

% Find ratio
rat = y_mlf1(3000)/y(3000) %#ok
figure;
plot(y_mlf1/rat); hold on;
plot(y);