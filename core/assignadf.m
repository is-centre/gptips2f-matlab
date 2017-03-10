function gp = assignadf(gp)
%ASSIGNADF Create ADF function handles in the caller workspace

% Get both function names and the corresponding evaluation symbols
funNames = gp.nodes.adf.name;
funSymbols = gp.nodes.adf.eval;

% Assign in caller's workspace
for k=1:length(funNames)
   assignin('caller', funNames{k}, str2func(funSymbols{k})); 
end

end

