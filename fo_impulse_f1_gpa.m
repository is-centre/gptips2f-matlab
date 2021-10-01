function y = fo_impulse_f1_gpa(t, a, alpha)
%FO_IMP_F1 Impulse response of a FO transfer function
%               1
%   G(s) = -----------
%          s^alpha + a

u = abs(a.*t.^alpha);
ppp = mlf_gpa_a(alpha,u);
figure;
subplot(211);
plot(u);
subplot(212);
plot(ppp);
y = t.^(alpha-1).*ppp;

% Different approach: just us a big value
% if isinf(y(1)), y(1)=1/eps; end

% Shift to the left, Right pad values
y = y(2:end); y(end+1)=y(end);

end

