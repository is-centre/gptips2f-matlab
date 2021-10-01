%% Define the original fotf
g = fotf('s^1.2+3s^0.4+5', 's^1.6+10s^1.2+35s^0.8+50s^0.4+24');

% Compute its impulse response
t = 0:0.001:2;
y1 = impulse(g,t);
yy1 = impulse(oustapp(g),t);

dt = t(2)-t(1);

%% Compute using MLF approximation
y2 = 71/6*t.^(-0.6).*mlf_a_a_app(0.4,4*t.^(0.4)) - ...
     31/2*t.^(-0.6).*mlf_a_a_app(0.4,3*t.^(0.4)) + ...
     9/2*t.^(-0.6).*mlf_a_a_app(0.4,2*t.^(0.4)) + ...
     1/6*t.^(-0.6).*mlf_a_a_app(0.4,t.^(0.4));
 
%% Step response
y2s = 71/6*t.^(0.4).*mlf_a_a1_app(0.4,4*t.^(0.4)) - ...
     31/2*t.^(0.4).*mlf_a_a1_app(0.4,3*t.^(0.4)) + ...
     9/2*t.^(0.4).*mlf_a_a1_app(0.4,2*t.^(0.4)) + ...
     1/6*t.^(0.4).*mlf_a_a1_app(0.4,t.^(0.4));
 
% Impulse?
y2i = colv(diff(y2s)/dt);
y2i = [y2i; y2i(end)];

%% Original MLF
y3 = 71/6*t.^(-0.6).*mlf(0.4,0.4,-4*t.^0.4) - ...
     31/2*t.^(-0.6).*mlf(0.4,0.4,-3*t.^0.4) + ...
     9/2*t.^(-0.6).*mlf(0.4,0.4,-2*t.^0.4) + ...
     1/6*t.^(-0.6).*mlf(0.4,0.4,-t.^0.4);

%% Plot it
figure; plot(t,y1); hold on;
plot(t,y2, 'r');
plot(t,y3, 'k');
plot(t,yy1);
plot(t,y2i);

%% figure;
plot(t,y2);

%% Try step response
global dt_conv;

u = ones(size(y1));
figure;
dt_conv = dt;
ys1 = convimp(u,y1);
ys2 = convimp(u,y2);
ys3 = convimp(u,y3);
ys4 = step(g,t);
ys5 = convimp(u,y2i);
plot(t,ys1); hold on;
plot(t,ys2,'x');
plot(t,ys3,'*');
plot(t,ys4);
plot(t,y2s);
plot(t,ys5);
xlim([0 2]);

legend('Orig impulse conv', 'MLF approx conv', ...
    'MLF conv', 'step in FOMCON', 'orig step mlfapprox', ...
    'step via mlfapprox conv imp...');


%% Now we try something completely different. I.e.: simulation
%  of the system response under arbitrary input signal.

global dt_conv;
dt_conv = 0.001;
st = refgen('prbs', 1, 50, 'Ts', dt_conv );
t = st.t;
plot(st);

%% Simulation via FOMCON
yf = lsim(g,st.u,st.t);

%% Simulation via step/impulse and convolution
ysi = 71/6*t.^(0.4).*mlf_a_a1_app(0.4,4*t.^(0.4)) - ...
     31/2*t.^(0.4).*mlf_a_a1_app(0.4,3*t.^(0.4)) + ...
     9/2*t.^(0.4).*mlf_a_a1_app(0.4,2*t.^(0.4)) + ...
     1/6*t.^(0.4).*mlf_a_a1_app(0.4,t.^(0.4));
 
% Impulse and step
yssi = colv(diff(ysi)/dt);
yssi = [yssi; yssi(end)];
qqq = conv(st.u, yssi);
yss = convimp(st.u, yssi);

%% Plot everything
figure; plot(st.t, yf); hold on; 
plot(st.t, yss);

%% Try a simple linear system
gg = tf(1,[1 1]);
yyyy = lsim(gg,st.u,st.t);

iiii = impulse(gg,st.t);
yyyy1 = 0.001*conv(st.u, iiii);
figure; plot(st.t, yyyy, st.t, yyyy1(1:length(st.t)));

