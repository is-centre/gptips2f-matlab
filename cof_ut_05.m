function y = cof_ut_05(u,t)
%COF_UT_05 Premade convolution-ready FO imp. response node

y = convimp(u,(fian_05(t)));

end

