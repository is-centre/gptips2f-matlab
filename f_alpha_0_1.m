function y = f_alpha_0_1(x)
%F_ALPHA_0_1 Returns the diff value of x^0.1*E_alpha,alpha+1(x^0.1)
%   (for impulse response)

%y = x.^0.9 .* mlf_gpa_a(0.1, abs(x).^0.1);
y = x.^0.1 .* mlf_a_a1_app(0.1, abs(x).^0.1);

% Only compute difference for 2+ elements
if numel(y) > 1
    y = diff(y);
    y(end+1) = y(end);
else
    % Do nothing
end

end

