% This script is used to generate examples
% commensurate-order FO systems with real rational poles

% Choose poles
pl = [-0.5 -3 -5 -9];

% Choose factors (for factorized form, which we will use)
fac = [1.5 10 1 -2];

% Choose the commensurate order (1 for classical systems)
alpha = 1;

% Generate the complete transfer function
full = 0;
for k=1:length(pl)
    % Show term and add it to the full expression
    term = fotf(num2str(fac(k)), ['s^' num2str(alpha) '+' num2str(-pl(k))]) %#ok
    full = full + term;
end

disp('Finally we get a system ...');
g = full %#ok

if (abs(alpha-1)<eps)
    g = oustapp(g);
end

t=0:0.01:10;
figure; step(g,t);