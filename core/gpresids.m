function varargout = gpresids(gp, ID, opts)
%VALIDATE Validate regression results
%
% Usage:  [ERR, STATS] = GPRESIDS(GP, ID, OPTS)
%
%         where ERR is the optional output argument which contains the
%                   absolute error vector y - y_reg
%               STATS is a structure with the following entries:
%                     .MaxError -- maximum absolute value of the simulation
%                                  error,
%                     .Mean -- mean value of residuals,
%                     .ResidNorm -- residual norm: sum(err^2)
%                     .MSE -- mean squared error: 1/N sum(err^2)
%                     .ConfBnd -- boundary value of the confidence interval
%                     .RErr -- vector with autocorrelation values of
%                              residuals for lags tau up until tau_max
%
%               If output arguments are omitted, plots validation results
%
%         Inputs:
%               GP is the GP data structure.
%               OPTS is a structure with the following entries (optional)
%                    .Knockout -- a boolean vector the with same number of
%                                 entries as genes in the individual to be
%                                 run. This evaluates the individual with
%                                 the indicated genes removed (in other
%                                 words, 'knocked out').
%                    .Conf -- measure of confidence in (0.0, 1.0]. Defaults
%                             to 0.95 (i.e., 95% confidence),
%                    .MaxTau -- maximum number of lags to consider for the
%                               autocorrelation of residuals test, defaults
%                               to 50

% Initialize output arguments
varargout = {};

% Number of arguments
if nargin < 2
    error('VALIDATE:NotEnoughInputArguments', ...
        'Not enough input arguments');
end

% Default values
maxtau = 50;
conf = 0.95;
knockout = 0;

% Check optional arguments
if nargin == 3
    
    if ~isa(opts, 'struct')
        error('RESIDS:OptionsNotAStructure', ...
            'The OPTS argument must be a structure');
    end
    
    % Knockout
    if gp_cfieldexists(opts, 'Knockout')
        knockout = opts.Knockout;
    end
    
    % Confidence
    if gp_cfieldexists(opts, 'Conf')
        conf = opts.Conf;
    end
    
    % Maximum lags
    if gp_cfieldexists(opts, 'MaxTau')
        maxtau = opts.MaxTau;
    end
    
end

if isempty(knockout) || ~any(knockout)
    doknockout = false;
else
    doknockout = true;
end

i = ID;

% Save graph state
showGraphState = gp.userdata.showgraphs;
gp.userdata.showgraphs = false; % Do not show graphs for this one

if isnumeric(ID)
    
    if i > 0 && i <= gp.runcontrol.pop_size
        
        %set this in case the fitness function needs to retrieve
        %the right returnvalues
        gp.state.current_individual = i;
        treestrs = tree2evalstr(gp.pop{i},gp);
        
        %if genes are being knocked out then remove appropriate gene
        if doknockout
            treestrs = kogene(treestrs, knockout);
            gp.state.force_compute_theta = true; %need to recompute gene weights if doing symbolic regression
        end
        
        % Here, use training data to compute the residuals
        residTxt = 'training data';
        y = gp.userdata.ytrain;
        [~,~,~,ypred] = feval(gp.fitness.fitfun,treestrs,gp);
        
    else
        error('A valid population member ID must be entered, e.g. 1, 99 or ''best''');
    end
    
elseif ischar(ID) && strcmpi(ID,'best')
    
    gp.fitness.returnvalues{gp.state.current_individual} = gp.results.best.returnvalues;
    treestrs = gp.results.best.eval_individual;
    
    if doknockout
        treestrs = kogene(treestrs, knockout);
        gp.state.force_compute_theta = true;
    end
    
    residTxt = 'training data';
    y = gp.userdata.ytrain;
    [~,~,~,ypred] = feval(gp.fitness.fitfun,treestrs,gp);
    
elseif ischar(ID) && strcmpi(ID,'valbest')
    
    % check that validation results/data present
    if ~isfield(gp.results,'valbest')
        disp('No validation results/data were found. Try runtree(gp,''best'') instead.');
        return;
    end
    
    %copy "valbest" return values to a slot in the "current" return values
    gp.fitness.returnvalues{gp.state.current_individual} = gp.results.valbest.returnvalues;
    
    treestrs = gp.results.valbest.eval_individual;
    
    if doknockout
        treestrs = kogene(treestrs, knockout);
        gp.state.force_compute_theta = true;
    end
    
    residTxt = 'validation data';
    y = gp.userdata.yval;
    [~,~,~,~,~,~,~,~,~,~,~,~,~,ypred]=feval(gp.fitness.fitfun,treestrs,gp);
    
elseif ischar(ID) && strcmpi(ID,'testbest')
    
    % check that test data results are present
    if ~isfield(gp.results,'testbest')
        disp('No test results/data were found. Try runtree(gp,''best'') instead.');
        return;
    end
    
    %copy "testbest" return values to a slot in the "current" return values
    gp.fitness.returnvalues{gp.state.current_individual} = gp.results.testbest.returnvalues;
    
    treestrs = gp.results.testbest.eval_individual;
    
    if doknockout
        treestrs = kogene(treestrs, knockout);
        gp.state.force_compute_theta = true;
    end
    
    residTxt = 'testing data';
    y = gp.userdata.ytest;
    [~,~,~,~,~,ypred]=feval(gp.fitness.fitfun,treestrs,gp);
    
    %if the selected individual is a GPMODEL struct (NB knockout disabled
    %for this form)
elseif isa(ID,'struct') && isfield(ID,'source') && ...
        (strcmpi(ID.source,'gpmodel2struct') || strcmpi(ID.source,'genes2gpmodel') );
    
    treestrs = ID.genes.geneStrs;
    treestrs = tree2evalstr(treestrs,gp);
    gp.fitness.returnvalues{gp.state.current_individual} = ID.genes.geneWeights;
    residTxt = 'training data';
    y = gp.userdata.ytrain;
    [~,~,~,ypred] = feval(gp.fitness.fitfun,treestrs,gp);
    
else
    error('Invalid argument.');
end

% maxtau must be less than the amount of simulated points
if length(y)<=maxtau
    maxtau = length(y)-1;
end

% Calculate simulation error
err = y - ypred;

% Confidence interval
intc = gp_quantile(c_p(conf))/sqrt(length(err));

% Compute necessary parameters
resnorm = sum(err.^2);              % Residual norm
mse     = 1/length(err)*resnorm;    % Mean squared error
errmean = mean(err);                % Mean value
[maxerror, ind] = max(abs(err));    % Maximum error

% Residuals at lags
rerr = [];
for k=0:maxtau
    rerr(end+1)=R_e(err,k);
end

% By first lag
rerr = rerr/rerr(1);

% Check no. of output arguments
if nargout == 0
    
    % Plot results
    h = figure();
    
    t = 1:length(err);
    
    subplot(2,1,1);
    plot(t, err, 'Color', 'r', 'Linewidth', 2);
    hold on;
    plot([t(1) t(end)], [errmean errmean], '--b', 'LineWidth', 2);
    stem(t(ind), err(ind), 'k', 'LineWidth', 2);
    xlabel('Sample');
    ylabel('Output error');
    title(['Mean squared error: ' num2str(mse) '; Max abs error: ' num2str(maxerror)]);
    grid;
    
    subplot(2,1,2);
    stem(1:maxtau,rerr(2:end), 'Color', 'b', 'Linewidth', 2);
    hold on;
    plot([1 maxtau], [intc intc], ':r', 'Linewidth', 2);
    plot([1 maxtau], -[intc intc], ':r', 'Linewidth', 2);
    xlabel('Lags [Samples]');
    xlim([1 maxtau]);
    title(['Autocorrelation of residuals (with P=' num2str(conf) ' confidence)']);
    
    set(h, 'NumberTitle', 'off');
    set(h, 'Name', ['Resid. analysis: ' residTxt]);
    
end

if nargout > 0
    varargout{1} = err;
end

if nargout > 1
    stats = struct;
    stats.MaxError = maxerror;
    stats.Mean = errmean;
    stats.ResidNorm = resnorm;
    stats.MSE = mse;
    stats.RErr = rerr(2:end); % Don't save the first one, it's value is 1.0
    stats.ConfBnd = intc;
    varargout{2} = stats;
end

% Restore graph state
gp.userdata.showgraphs = showGraphState;

end

% Residual covariance
function Re = R_e(err,tau)
Re = 1/(length(err)-tau)*sum(err(1:end-tau).*err(1+tau:end));
end

% Argument for quantile function: assuming normal distribution
function cp = c_p(p)
cp = 1-0.5*(1-p);
end

function y = gp_quantile(x)
%QUANTILE Compute quantile function Phi^-1(x)
    y = sqrt(2)*erfinv(2*x-1);
end


function [fexist, f] = gp_cfieldexists(structure, structure_fields)
%CFIELDEXISTS Check if a field (given by a string) in a structure exists

fexist = 1;  % Unless proven otherwise
fields_to_check = explode(structure_fields, '.');
k=1;

substructure = structure;
while (fexist && k<=length(fields_to_check))
    % Not a structure or Field not found
    if ~isstruct(substructure) || ...
            ~(any(strcmp(fieldnames(substructure), fields_to_check{k})))
        fexist = 0;
    else
        substructure = substructure.(fields_to_check{k});
    end
    k=k+1;
end

% Return the field, if requested
if fexist
    f = substructure;
else
    f = [];
end

end

