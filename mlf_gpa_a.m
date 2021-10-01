function y = mlf_gpa_a(alpha, x)
%MLF_GPA_A Computes Mittag-Leffler function global pade approximation
%
% Here, we consider the case E     (-x), i.e. single parameter (alpha) MLF
%                             alpha

% Use the Pade approximation
g_1ma = gamma(1-alpha);
g_1pa = gamma(1+alpha);
g_1m2a = gamma(1-2*alpha);

% Elements of the polynomials
bx = [1/gamma(alpha)];
ax = [g_1ma/g_1pa 2*g_1ma^2/(g_1pa*g_1m2a) 1];

% Construct the result
y = polyval(bx,x)./polyval(ax,x);

end

