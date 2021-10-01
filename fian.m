function y = fian(x,alpha)
%FIMP_A_N Returns part of impulse response of a FO system 1/(s^alpha+a)

% Check second input: must be in valid range. This depends on the
% implementation. Currently we allow for .1 resolution in the range .1-.9.

% There's a whole bunch of operations going on here.
% Make sure that you interpret them correctly in the output function!
alpha = floor(10*abs(alpha(1)))/10;
if (alpha < .1) alpha = .1; end % TODO: Should do this check centrally
if (alpha > .9) alpha = .9; end %  through some related function

% We consider the absolute value of x, otherwise results will be erroneous
x = abs(x);

% We can monitor alphas! (be careful to reset this regularly...)
% To access, call "global alphas" from Workspace
% global alphas;
% alphas(end+1) = alpha;

global dt_conv;

% Prolong x by a single element
if numel(x) > 1
    dt = x(end)-x(end-1);
else
    dt = x(end);
end

% Automatically compensate for n when the input argument is n*t
if ~(abs(dt - dt_conv)<eps)
    cfac = 1/((dt/dt_conv)^alpha);
else
    cfac = 1;
end

x(end+1) = x(end) + dt;

y = cfac.* x.^alpha .* mlf_a_a1_app(alpha, x.^alpha);
y = diff(y)/dt_conv;

end

