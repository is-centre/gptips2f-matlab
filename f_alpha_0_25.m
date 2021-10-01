function y = f_alpha_0_25(x)
%F_ALPHA_0_1 Returns the value of x^0.75*E_alpha,alpha_approx(x^0.25)
%   E_alpha,alpha_approx() is the approximated Mittag-Leffler function

y = x.^0.75 .* mlf_gpa_a(0.25, abs(x).^0.25);

end

