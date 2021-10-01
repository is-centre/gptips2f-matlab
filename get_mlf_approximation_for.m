% This script computes a Pade approximation for the specified MLF

% Let's get started with the following: MLF(0.5,0.5,x), where x spans
% from 10e-5 to 10e5; the number of points can be modified for better
% precision. We try N=100000.

% Parameters
alpha = 0.9;
x = linspace(10^(-5),10^5,10^7);
prec = 10;

% Get the response from the original function
tim = tic;
y = mlf(alpha,alpha,-x,prec);
toc(tim)

% See the result
figure; plot(x,y);

%% Now the hard part is over ... optimization time!

%% Set the polynomial order
zer = 2;
pol = 5;

% And corresponding ranges
rzer = 1:zer;
rpol = (1+zer):(zer+pol);

%% The rational polynomial is comprised of zer+pol parameters
p = ones(zer+pol,1);

%% Start from where we left?
p = xo;

%% Do the iterations...

% Let us define the objective: sum of squares
% f = @(z) sum((1./x) .* (y - polyval(z(1:5),-x)./polyval(z(6:10),-x)).^2);

% Or ISE
f = @(z) trapz(x, (y - polyval(z(rzer),-x)./polyval(z(rpol),-x)).^2);

% For lsqnonlin
f1 = @(z) y - polyval(z(rzer),-x)./polyval(z(rpol),-x);

% This one to test the actual values
f_ = @(z,c) polyval(z(rzer),-c) ./ polyval(z(rpol), -c);

%% We'll be using our good friend fminsearch() initially
xo = fminsearch(f, p, optimset('Display', 'iter', 'MaxFunEvals', 20000, ...
    'MaxIter', 10000));

%% And eventually we'll try nonlinear least squares as well...
xo = lsqnonlin(f1, p, [], [], optimset('Display', 'iter', 'MaxFunEvals', 20000, ...
    'MaxIter', 10000, 'Algorithm', 'levenberg-marquardt'));

%% Evaluate result: show error

% Compute the abs err sum as well
err = f_(xo,x)-y;

% This is measurement of quality against original MLF
err_sum = sum(abs(err)); 

figure; 
h1 = subplot(211);
plot(x,y); hold on; plot(x,f_(xo,x));
h2 = subplot(212);
plot(x,f_(xo,x)-y); title(['Approximation error; absolute err sum: ' ...
    num2str(err_sum)]);
linkaxes([h1, h2], 'x');

%% Evaluate result with oo: show error
figure; 
h1 = subplot(211);
plot(x,y); hold on; plot(x,f_(oo,x));
h2 = subplot(212);
plot(x,f_(oo,x)-y); title('Approximation error');
linkaxes([h1, h2], 'x');