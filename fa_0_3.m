function y = fa_0_3(x)
%FA_0_3 Returns the imp response of a system with alpha=0.3

global dt_conv; % This is still necessary!

% Prolong x by a single element
if numel(x) > 1
    dt = x(end)-x(end-1);
else
    dt = x(end);
end

x(end+1) = x(end) + dt;

y = x.^0.3 .* mlf_a_a1_app(0.3, abs(x).^0.3);
y = diff(y)/dt_conv;

end

