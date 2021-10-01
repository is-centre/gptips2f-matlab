function y = f_alpha_0_4(x)
%F_ALPHA_0_4 Returns the diff value of x^0.4*E_alpha,alpha_approx(x^0.4)

global dt_conv;

y = x.^0.4 .* mlf_a_a1_app(0.4, abs(x).^0.4);

% Only compute difference for 2+ elements
if numel(y) > 1
    y = diff(y)/dt_conv;
    y(end+1) = y(end);
else
    % Do nothing
end

end

