%% [TEST OK] Integer-order system example
t=0:0.001:2; u=ones(size(t)); % dt appears to be of critical importance
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['71/6*convimp(u,exp(-4*t))' ...
        '-31/2*convimp(u,exp(-3*t))' ...
        '+9/2*convimp(u,exp(-2*t))' ...
        '+1/6*convimp(u,exp(-1*t))' ...
        ]);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y,'r');

% Compare to "correct" response
g = fotf('s^3+3s+5', 's^4+10s^3+35s^2+50s+24');
hold on; step(oustapp(g),t);

%% [TEST OK] Alternate integer-order system example
t=0:0.001:2; u=ones(size(t)); % dt appears to be of critical importance
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['71/6*co_ut(u,-4*t)' ...
        '-31/2*co_ut(u,-3*t)' ...
        '+9/2*co_ut(u,-2*t)' ...
        '+1/6*co_ut(u,-1*t)' ...
        ]);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y,'r');

% Compare to "correct" response
g = fotf('s^3+3s+5', 's^4+10s^3+35s^2+50s+24');
hold on; step(oustapp(g),t);


%% [TEST OK] Equivalence test with constants
t=0:0.01:2; u=ones(size(t));
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['71/6*convimp(u,(fian_04(4^(1/0.4)*t)))' ...
        '-31/2*convimp(u,(fian_04(3^(1/0.4)*t)))' ...
        '+9/2*convimp(u,(fian_04(2^(1/0.4)*t)))' ...
        '+1/6*convimp(u,(fian_04(1^(1/0.4)*t)))' ...
        ]);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y)

%% [TEST OK] Repeated poles example: GENERALIZED VERSION
t=0:0.01:2; u=ones(size(t));
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['-convimp(u,(fian(2^(1/0.3)*t, 0.3)))' ...
        '+convimp(u,(fian(1^(1/0.3)*t, 0.3)))' ...
        '-convimp(u,convself((fian(1^(1/0.3)*t, 0.3))))' ...
        '+convimp(u,convself2((fian(1^(1/0.3)*t, 0.3))))']);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y)

%% No repeated pole example
t=0:0.01:2; u=ones(size(t));
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['convimp(u,71/6*(fian_04(4^(1/0.4)*t)))' ...
        '-convimp(u,31/2*(fian_04(3^(1/0.4)*t)))' ...
        '+convimp(u,9/2*(fian_04(2^(1/0.4)*t)))' ...
        '+convimp(u,1/6*(fian_04(1^(1/0.4)*t)))' ...
        ]);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y)


%% Repeated poles example
t=0:0.01:2; u=ones(size(t));
global dt_conv;
dt_conv = t(2)-t(1);

solf = sym(['-convimp(u,1/2*(fa_0_3(2^(1/0.3)*t)))' ...
        '+convimp(u,(fa_0_3(1^(1/0.3)*t)))' ...
        '-convimp(u,convself((fa_0_3(1^(1/0.3)*t))))' ...
        '+convimp(u,convself2((fa_0_3(1^(1/0.3)*t))))']);
    
myTest = matlabFunction(solf);

y=myTest(t,u);
figure; plot(t,y)

