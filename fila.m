function a = fila(a)
%FILA Filter the value of A such that it is inside the range [0.1, 0.9]

% Must be of proper format
a = str2double(sprintf('%.1f',a));

% Saturation
if (a<0.1||isnan(a)) a=0.1; end;
if (a>0.9) a=0.9; end;

end

