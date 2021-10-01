% The idea of this script is to test the residue
% approach, i.e. partial fraction expansion of the
% pseudo rational transfer function of a FOTF. With
% this approach, we can easily use the MLF approximations
% to compute step responses of dynamic systems. For now,
% we limit the class of the studied systems to those without
% repeated and/or complex poles in the pseudo-rational pole polynomial.


%% Define the original fotf
g = fotf('s^1.2+3s^0.4+5', 's^1.6+10s^1.2+35s^0.8+50s^0.4+24')

%% Another one
g = fotf('1', 's^1.2+5s^0.9+9s^0.6+7s^0.3+2')

%% One more
g = fotf('1','s^{1.8}+s^{0.9}+1');

%% This one will be interesting
loadsets;
g = G2

%% Do the stuff
t=0:0.001:100; % Time vector

[b,a,alpha] = tfdata(g)

% Try residue
[r,p,k] = residue(b,a)

% Compute step response assuming NO REPEATED or COMPLEX poles
if numel(unique(p))<numel(p)
    error('Cannot handle repeated poles.');
end

if any(imag(r)) || any(imag(p))
    error('Cannot handle complex zeros or poles.');
end


resp = zeros(length(t),1);
for k=1:numel(p)
    smd = t.^alpha.*r(k).*mlf_a_a1_app(alpha,-p(k).*t.^alpha);
    resp = resp + colv(smd);
end

figure; plot(t,resp);

%% The alternate test
global dt_conv;
dt_conv = t(2)-t(1);
us = ones(size(t));
resp = zeros(length(t),1);
for k=1:numel(p)
    smd = 1/abs(p(k)).*r(k).*f_alpha_0_4(abs(p(k))^(1/0.4).*t);
    resp = resp + convimp(us, colv(smd));
end
hold on; plot(t,resp);

%% The alternate test V2: response to arbitrary input
st = refgen('prbs', 1, 50, 'Ts', dt_conv );
t = st.t;
us = st.u;
global dt_conv;
dt_conv = t(2)-t(1);
% us = ones(size(t));
resp = zeros(length(t),1);
for k=1:numel(p)
    smd = 1/abs(p(k)).*r(k).*f_alpha_0_4(abs(p(k))^(1/0.4).*t);
    resp = resp + convimp(us, colv(smd));
end
hold on; plot(t,resp);

%% Test for a system with repeated poles
Gg = fotf('1','(s^0.5+2)^3');
Go = fotf('1','s^0.5+2');

% Technically, to get Go to behave like Gg, we need
% to convolve the response with itself. Should work
% for both step and impulse responses.
t=0:0.01:100;
dt = t(2)-t(1);
global dt_conv;
dt_conv = dt;
y1 = step(Gg,t);
y1 = y1(1:length(t));
figure; plot(t,y1);

%% What about our non-pole-repeated friend?
y21 = step(Go,t);
hold on; plot(t,y21); % Looks a tad similar? Let's convolve it with itself
y22 = dt*conv(y21,y21);

%% How does it look?
plot(y22); % Like shit. Let's try impulse response approach.

%%
y31 = impulse(Go,t); % Now we'll convolve THIS with itself and check it out
figure; plot(t,y31);

%% Another way of impulse ... ^3 now
y31 = t.^(0.5).*mlf_a_a1_app(0.5,2*t.^(0.5));
y31 = diff(y31)./dt;
y31 = colv(y31, y31(end));
hold on; plot(t,y31);

%% Convolve as many times as the pole is repeated
y33 = convself2(y31);

%% Now, convolve with ones
y3 = dt*conv(ones(length(t),1), y33);
y3 = y3(1:length(t));
figure; plot(y1); hold on; plot(y3);
