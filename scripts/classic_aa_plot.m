function [matSamLat, matLatSam, obj] = classic_aa_plot(matFeatSam, nLat)
% This function calls an R script to run standard archetypal analysis
%
% Supporting function: classic_aa_plot.R
% Intermediate files: classic_aa_Input.mat classic_aa_Output.mat (deleted after run by default)
%
% Note: change path to Rscript accordingly
%
% Copyright Sohan Seth sohan.seth@hiit.fi

save classic_aa_Input.mat matFeatSam nLat
% In triton
% [status, msg] = system('/m/fs/software/R/ubuntu/R-2.15.1/bin/Rscript classic_aa.R');
% [status, msg] = system('Rscript classic_aa.R');
% if status ~= 0
%     disp(msg);
% end
!/opt/local/bin/Rscript classic_aa_plot.R
load classic_aa_Output.mat matSamLat matLatSam obj
!rm classic_aa_Input.mat classic_aa_Output.mat