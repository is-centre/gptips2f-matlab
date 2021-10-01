% This one takes care of repeated poles problem ...

%% Another one
g = fotf('1', 's^1.2+5s^0.9+9s^0.6+7s^0.3+2')

%% Do the stuff
t=0:0.001:2; % Time vector
global dt_conv;
dt_conv = t(2)-t(1);

[b,a,alpha] = tfdata(g)

% Try residue
[r,p,k] = residue(b,a)

if any(imag(r)) || any(imag(p))
    error('Cannot handle complex zeros or poles.');
end

% Impulse response
iresp = zeros(length(t),1);

% We'll need to construct a map of the poles to check for repeated poles
pmap = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k=1:length(p)
    
    % First we compute the step response as usual for a single component
    smd = t.^alpha.*mlf_a_a1_app(alpha,abs(p(k)).*t.^alpha);

    % Differentiate to get impulse response
    dsmd = diff(smd)./dt_conv; dsmd = colv(dsmd, dsmd(end)); % Repeat last value
    
    % We now check if the signal needs to be convolved with itself
    if (~isKey(pmap, p(k)))
        pmap(p(k)) = 1;
    else
        dsmd = convselfn(dsmd, pmap(p(k)));
        pmap(p(k)) = pmap(p(k))+1;
    end
 
    iresp = iresp + r(k).*dsmd;
    
end

% Finally, just convolve with 1-s to get step response
u = ones(length(t),1);
y = convimp(u, iresp);
y1 = step(g,t);

% Error
figure; plot(t,y); hold on; plot(t,);

