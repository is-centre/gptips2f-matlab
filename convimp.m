function y = convimp(u, v)
%CONVIMP Convolve two uniformly sampled signal vectors
%   A global variable is dt_conv is needed for correct computation
%   such that holds the fixed step, i.e., sampling interval, in seconds.

% Default, fallback value or dt
DEFAULT_DT = 1e-3;

% Convert to column vectors if needed
u = colv(u);
v = colv(v);

global dt_conv
if (~exist('dt_conv', 'var') || isempty(dt_conv))
    dt_conv = DEFAULT_DT;
end

y = conv(u,v);
y = dt_conv * y(1:length(u));

% This is a test feature that ensures y(1) == 0
% TODO: Note that here we assume zero initial conditions and t(1) == 0
y(1) = 0;

end

