function y = cof_ut_01(u,t)
%COF_UT_01 Premade convolution-ready FO imp. response node

y = convimp(u,(fian_01(t)));

end

