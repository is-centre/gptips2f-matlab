function y = convselfnt(u,n)
%CONVSELFN Convolve signal with itself N times

% GPTIPS version. Should be careful!
n = floor(abs(n(1))); % Always use first element regardless of dimension

% Because of the floor operation, the corresponding 
% expression should be interpreted accordingly!

% Default, fallback value or dt
DEFAULT_DT = 1e-3;

global dt_conv
if (~exist('dt_conv', 'var') || isempty(dt_conv))
    dt_conv = DEFAULT_DT;
end

% Convert to column vector
u = colv(u);

if n==0
    y = u;
elseif n>0
    for k=1:n
        if k==1
            y = dt_conv * conv(u,u);
        else
            y = dt_conv * conv(u,y);
        end
        y =  y(1:length(u));
    end
end

end


