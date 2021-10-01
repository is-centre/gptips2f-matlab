g = fotf('1','s^0.3');

expr = @(t) 3.3333*t.^(0.3)./gamma(0.3);

% Try the digital thing
sum = 0;

%%
for k=1:t
    
end

%%
t=1e-5:0.01:1000;


%%
step(g,t);
hold on;
plot(t,expr(t));