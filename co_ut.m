function y = co_ut(u,t)
%CO_UT Premade convolution-ready FO imp. response node

y = convimp(u,(exp(t)));

end

