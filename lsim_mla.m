function [y, yi] = lsim_mla(G,u,t)
%LSIM_MLA Linear simulation of commensurate-order FOTF systems
%   [Y, YI] = LSIM_MLA(G,U,T) computes the response of a FOTF system in G,
%                             under input signal U, at time instances in T.
%                             The impulse response YI is also returned.

% Get fotf data and commensurate order
[b,a,alpha] = tfdata(G);

% Partial fraction expansion
[r,p,k] = residue(b,a);

% TODO: direct terms
if ~isempty(k)
    error('Direct terms resulting in PFE not supported.');
end

% TODO: add some support for complex poles
if any(imag(r)) || any(imag(p))
    error('Cannot handle complex zeros or poles.');
end

% Impulse response
iresp = zeros(length(t),1);

% Assuming regular sampling interval
dt = t(2)-t(1);

% We'll need to construct a map of the poles to check for repeated poles
pmap = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k=1:length(p)
    
    % First we compute the step response as usual for a single component
    % We add an extra entry for t to accomodate for the impulse response
    t1 = colv(t, t(end)+dt);
    smd = t1.^alpha.*mlf_a_a1_app(alpha,abs(p(k)).*t1.^alpha);

    % Differentiate to get impulse response
    dsmd = diff(smd)./dt;
    
    % We now check if the signal needs to be convolved with itself
    if (~isKey(pmap, p(k)))
        pmap(p(k)) = 1;
    else
        r(k)
        dsmd = convselfn(dsmd, pmap(p(k)), dt);
        pmap(p(k)) = pmap(p(k))+1;
    end
 
    % Add component of impulse response
    iresp = iresp + r(k).*dsmd;
    
end

% Finally, just convolve with u(t) to get response
y = convimp(u, iresp);
yi = iresp;

end

