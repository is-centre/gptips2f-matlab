function E = ml(z,alpha,beta,gama)
%
% Evaluation of the Mittag-Leffler (ML) function with 1, 2 or 3 parameters
% by means of the OPC algorithm [1]. The routine evaluates an approximation
% Et of the ML function E such that |E-Et|/(1+|E|) approx 1.0e-15   
%     
%
% E = ML(z,alpha) evaluates the ML function with one parameter alpha for
% the corresponding elements of z; alpha must be a real and positive
% scalar. The one parameter ML function is defined as
%
% E = sum_{k=0}^{infty} z^k/Gamma(alpha*k+1)
%
% with Gamma the Euler's gamma function.
%
%
% E = ML(z,alpha,beta) evaluates the ML function with two parameters alpha
% and beta for the corresponding elements of z; alpha must be a real and
% positive scalar and beta a real scalar. The two parameters ML function is
% defined as
%
% E = sum_{k=0}^{infty} z^k/Gamma(alpha*k+beta)
%
%
% E = ML(z,alpha,beta,gama) evaluates the ML function with three parameters
% alpha, beta and gama for the corresponding elements of z; alpha must be a
% real scalar such that 0<alpha<1, beta any real scalar and gama a real and
% positive scalar; the arguments z must satisfy |Arg(z)| > alpha*pi. The
% three parameters ML function is defined as 
%
% E = sum_{k=0}^{infty} Gamma(gama+k)*z^k/Gamma(gama)/k!/Gamma(alpha*k+beta)
%
%
% NOTE: 
% This routine implements the optimal parabolic contour (OPC) algorithm
% described in [1] and based on the inversion of the Laplace transform on a
% parabolic contour suitably choosen in one of the regions of analyticity
% of the Laplace transform.
%
%
% REFERENCES
%
%   [1] R. Garrappa, Numerical evaluation of two and three parameter
%   Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 2015,
%   53(3), 1350-1369
%
%
%   Please, report any problem or comment to : 
%        roberto dot garrappa at uniba dot it
%
%   Copyright (c) 2015, Roberto Garrappa, University of Bari, Italy
%   roberto dot garrappa at uniba dot it
%   Homepage: http://www.dm.uniba.it/Members/garrappa
%   Revision: 1.4 - Date: October 8 2015


% Check inputs
if nargin < 4
    gama = 1 ;
    if nargin < 3
        beta = 1 ;
        if nargin < 2
            error('MATLAB:ml:NumberParameters', ...
                'The parameter ALPHA must be specified.');
        end
    end
end

% Check whether the parameters ALPHA, BETA and GAMMA are scalars
if length(alpha) > 1 || length(beta) > 1 || length(gama) > 1
    alpha = alpha(1) ; beta = beta(1) ; gama = gama(1) ;
    warning('MATLAB:ml:ScalarParameters', ...
        ['ALPHA, BETA and GAMA must be scalar parameters. ', ...
        'Only the first values ALPHA=%f BETA=%f and GAMA=%f will be used. '], ...
        alpha, beta, gama) ;
end

% Check whether the parameters meet the contraints
if real(alpha) <= 0 || real(gama) <= 0 || ~isreal(alpha) || ...
        ~isreal(beta) || ~isreal(gama)
    error('MATLAB:ml:ParametersOutOfRange', ...
        ['Error in the parameters of the Mittag-Leffler function. ', ...
        'Parameters ALPHA and GAMA must be real and positive. ', ...
        'The parameter BETA must be real.']) ;
end

% Check parameters and arguments for the three parameter case
if abs(gama-1) > eps
    if alpha > 1
        error('MATLAB:ml:ALPHAOutOfRange',...
            ['With the three parameters Mittag-Leffler function ', ...
            'the parameter ALPHA must satisfy 0 < ALPHA < 1']) ;
    end
    if min(abs(angle(z(abs(z)>eps)))) <= alpha*pi 
        error('MATLAB:ml:ThreeParametersArgument', ...
            ['With the three parameters Mittag-Leffler function ', ...
            'this code works only when |Arg(z)|>alpha*pi.']);
    end
end

% Target precision 
log_epsilon = log(10^(-15)) ; 

% Inversion of the LT for each element of z
E = zeros(size(z)) ;  
for k = 1 : length(z)
    if abs(z(k)) < 1.0e-15
        E(k) = 1/gamma(beta) ; 
    else
        E(k) = LTInversion(1,z(k),alpha,beta,gama,log_epsilon) ;
    end
end


end


% =========================================================================
% Evaluation of the ML function by Laplace transform inversion
% =========================================================================
function E = LTInversion(t,lambda,alpha,beta,gama,log_epsilon)

% Evaluation of the relevant poles
theta = angle(lambda) ;
kmin = ceil(-alpha/2 - theta/2/pi) ;
kmax = floor(alpha/2 - theta/2/pi) ;
k_vett = kmin : kmax ;
s_star = abs(lambda)^(1/alpha) * exp(1i*(theta+2*k_vett*pi)/alpha) ;

% Evaluation of phi(s_star) for each pole
phi_s_star = (real(s_star)+abs(s_star))/2 ;

% Sorting of the poles according to the value of phi(s_star)
[phi_s_star , index_s_star ] = sort(phi_s_star) ;
s_star = s_star(index_s_star) ;

% Deleting possible poles with phi_s_star=0
index_save = phi_s_star > 1.0e-15 ;
s_star = s_star(index_save) ;
phi_s_star = phi_s_star(index_save) ;

% Inserting the origin in the set of the singularities
s_star = [0, s_star] ;
phi_s_star = [0, phi_s_star] ;
J1 = length(s_star) ; J = J1 - 1 ;

% Strength of the singularities
p = [ max(0,-2*(alpha*gama-beta+1)) , ones(1,J)*gama ]  ;
q = [ ones(1,J)*gama , +Inf] ;
phi_s_star = [phi_s_star, +Inf] ;

% Looking for the admissible regions with respect to round-off errors
admissible_regions = find( ...
    (phi_s_star(1:end-1) < (log_epsilon - log(eps))/t) & ...
    (phi_s_star(1:end-1) < phi_s_star(2:end))) ;

% Initializing vectors for optimal parameters
JJ1 = admissible_regions(end) ;
mu_vett = ones(1,JJ1)*Inf ;
N_vett = ones(1,JJ1)*Inf ;
h_vett = ones(1,JJ1)*Inf ;

% Evaluation of parameters for inversion of LT in each admissible region
find_region = 0 ;
while ~find_region
    for j1 = admissible_regions
        if j1 < J1
            [muj,hj,Nj] = OptimalParam_RB ...
                (t,phi_s_star(j1),phi_s_star(j1+1),p(j1),q(j1),log_epsilon) ;
        else
            [muj,hj,Nj] = OptimalParam_RU(t,phi_s_star(j1),p(j1),log_epsilon) ;
        end
        mu_vett(j1) = muj ; h_vett(j1) = hj ; N_vett(j1) = Nj ;
    end
    if min(N_vett) > 200
        log_epsilon = log_epsilon +log(10) ;
    else
        find_region = 1 ;
    end
end
    

% Selection of the admissible region for integration which involves the
% minimum number of nodes 
[N, iN] = min(N_vett) ; mu = mu_vett(iN) ; h = h_vett(iN) ;

% Evaluation of the inverse Laplace transform
k = -N : N ;
u = h*k ;
z = mu*(1i*u+1).^2 ;
zd = -2*mu*u + 2*mu*1i ;
zexp = exp(z*t) ;
F = z.^(alpha*gama-beta)./(z.^alpha - lambda).^gama.*zd ;
S = zexp.*F ;
Integral = h*sum(S)/2/pi/1i ;

% Evaluation of residues
ss_star = s_star(iN+1:end) ;
Residues = sum(1/alpha*(ss_star).^(1-beta).*exp(t*ss_star)) ;

% Evaluation of the ML function
E = Integral + Residues ;
if isreal(lambda) 
    E = real(E) ;
end

end


% =========================================================================
% Finding optimal parameters in a right-bounded region
% =========================================================================
function [muj,hj,Nj] = OptimalParam_RB ...
    (t, phi_s_star_j, phi_s_star_j1, pj, qj, log_epsilon)

% Definition of some constants
log_eps = -36.043653389117154 ; % log(eps)
fac = 1.01 ;
conservative_error_analysis = 0 ;

% Maximum value of fbar as the ration between tolerance and round-off unit
f_max = exp(log_epsilon - log_eps) ;

% Evaluation of the starting values for sq_phi_star_j and sq_phi_star_j1
sq_phi_star_j = sqrt(phi_s_star_j) ;
threshold = 2*sqrt((log_epsilon - log_eps)/t) ;
sq_phi_star_j1 = min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j) ;

% Zero or negative values of pj and qj
if pj < 1.0e-14 && qj < 1.0e-14
    sq_phibar_star_j = sq_phi_star_j ;
    sq_phibar_star_j1 = sq_phi_star_j1 ;
    adm_region = 1 ;
end

% Zero or negative values of just pj
if pj < 1.0e-14 && qj >= 1.0e-14
    sq_phibar_star_j = sq_phi_star_j ;
    if sq_phi_star_j > 0
        f_min = fac*(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j))^qj ;
    else
        f_min = fac ;
    end
    if f_min < f_max
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fq = f_bar^(-1/qj) ;
        sq_phibar_star_j1 = (2*sq_phi_star_j1-fq*sq_phi_star_j)/(2+fq) ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end

% Zero or negative values of just qj
if pj >= 1.0e-14 && qj < 1.0e-14
    sq_phibar_star_j1 = sq_phi_star_j1 ;
    f_min = fac*(sq_phi_star_j1/(sq_phi_star_j1-sq_phi_star_j))^pj ;
    if f_min < f_max
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fp = f_bar^(-1/pj) ;
        sq_phibar_star_j = (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp) ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end

% Positive values of both pj and qj
if pj >= 1.0e-14 && qj >= 1.0e-14
    f_min = fac*(sq_phi_star_j+sq_phi_star_j1)/...
        (sq_phi_star_j1-sq_phi_star_j)^max(pj,qj) ;
    if f_min < f_max
        f_min = max(f_min,1.5) ;
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fp = f_bar^(-1/pj) ;
        fq = f_bar^(-1/qj) ;
        if ~conservative_error_analysis
            w = -phi_s_star_j1*t/log_epsilon ;
        else
            w = -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t) ;
        end
        den = 2+w - (1+w)*fp + fq ;
        sq_phibar_star_j = ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den ;
        sq_phibar_star_j1 = (-(1+w)*fq*sq_phi_star_j ...
            + (2+w-(1+w)*fp)*sq_phi_star_j1)/den ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end

if adm_region
    log_epsilon = log_epsilon  - log(f_bar) ;
    if ~conservative_error_analysis
        w = -sq_phibar_star_j1^2*t/log_epsilon ;
    else
        w = -2*sq_phibar_star_j1^2*t/(log_epsilon-sq_phibar_star_j1^2*t) ;
    end
    muj = (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2 ;
    hj = -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)...
        /((1+w)*sq_phibar_star_j + sq_phibar_star_j1) ;
    Nj = ceil(sqrt(1-log_epsilon/t/muj)/hj) ;
else
    muj = 0 ; hj = 0 ; Nj = +Inf ;
end

end

% =========================================================================
% Finding optimal parameters in a right-unbounded region
% =========================================================================
function [muj,hj,Nj] = OptimalParam_RU (t, phi_s_star_j, pj, log_epsilon)

% Evaluation of the starting values for sq_phi_star_j
sq_phi_s_star_j = sqrt(phi_s_star_j) ;
if phi_s_star_j > 0
    phibar_star_j = phi_s_star_j*1.01 ;
else
    phibar_star_j = 0.01 ;
end
sq_phibar_star_j = sqrt(phibar_star_j) ;

% Definition of some constants
f_min = 1 ; f_max = 10 ; f_tar = 5 ;

% Iterative process to look for fbar in [f_min,f_max]
stop = 0 ;
while ~stop
    phi_t = phibar_star_j*t ; log_eps_phi_t = log_epsilon/phi_t ;
    Nj = ceil(phi_t/pi*(1 - 3*log_eps_phi_t/2 + sqrt(1-2*log_eps_phi_t))) ;
    A = pi*Nj/phi_t ;
    sq_muj = sq_phibar_star_j*abs(4-A)/abs(7-sqrt(1+12*A)) ;
    fbar = ((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)^(-pj) ;
    stop = (pj < 1.0e-14) || (f_min < fbar && fbar < f_max) ;
    if ~stop
        sq_phibar_star_j = f_tar^(-1/pj)*sq_muj + sq_phi_s_star_j ;
        phibar_star_j = sq_phibar_star_j^2 ;
    end
end
muj = sq_muj^2 ;
hj = (-3*A - 2 + 2*sqrt(1+12*A))/(4-A)/Nj ;

% Adjusting integration parameters to keep round-off errors under control
log_eps = log(eps) ; threshold = (log_epsilon - log_eps)/t ;
if muj > threshold
    if abs(pj) < 1.0e-14 , Q = 0 ; else Q = f_tar^(-1/pj)*sqrt(muj) ; end
    phibar_star_j = (Q + sqrt(phi_s_star_j))^2 ;
    if phibar_star_j < threshold
        w = sqrt(log_eps/(log_eps-log_epsilon)) ;
        u = sqrt(-phibar_star_j*t/log_eps) ;
        muj = threshold ;
        Nj = ceil(w*log_epsilon/2/pi/(u*w-1)) ;
        hj = sqrt(log_eps/(log_eps - log_epsilon))/Nj ;
    else
        Nj = +Inf ; hj = 0 ;
    end
end

end

