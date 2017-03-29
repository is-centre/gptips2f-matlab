function toolboxes = gp_toolboxlist()
%GP_TOOLBOXLIST Generate a list of currently available MATLAB toolboxes

v = ver;
[toolboxes{1:length(v)}] = deal(v.Name);

end

