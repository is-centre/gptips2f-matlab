function y = cof_ut_04(u,t)
%COF_UT_04 Premade convolution-ready FO imp. response node

y = convimp(u,(fian_04(t)));

end

