function y = cof_ut_03(u,t)
%COF_UT_03 Premade convolution-ready FO imp. response node

y = convimp(u,(fian_03(t)));

end

