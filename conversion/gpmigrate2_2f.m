function gp = gpmigrate2_2f(gp)
%GPMIGRATE2_2F Migrate parameters of GP structure from V2 to V2F of toolbox
%   Some parameters of the GP strucure have changed in the new version of
%   the toolbox. Therefore, this function is provided to update those
%   parameters. Simply do GP = GPMIGRATE2_2F(GP) to get those parameters
%   fixed.
%   
%   List of parameter changes:
%       * Toolboxes are now correctly detected on all computers
%       dynamically by means of anonymous function handles. This means that
%       when you run GPTIPS 2F on a computer that doesn't have MATLAB's
%       Symbolic Math toolbox installed, for instance, and then use the
%       resulting GP structure on a computer that does have it installed,
%       the toolbox is detected and you can proceed.

% Correct toolbox check using handles
[gp.info.toolbox.symbolic, gp.info.toolbox.parallel, gp.info.toolbox.stats] = gptoolboxcheck;

end

