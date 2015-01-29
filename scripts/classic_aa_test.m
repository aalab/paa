function matLatSam = classic_aa_test(matFeatSam)
% This function calls an R script to estimate H using existing archetypes
%
% Supporting function: classic_aa_test.R
% Intermediate files: classic_aa_test_Input.mat classic_aa_test_Output.mat (deleted after run by default)
% Archetypes are assumed to exist in bestModel.Rdat
%
% Note: change path to Rscript accordingly
%
% Copyright Sohan Seth sohan.seth@hiit.fi

save classic_aa_test_Input.mat matFeatSam
% In triton
% [status, msg] = system('/m/fs/software/R/ubuntu/R-2.15.1/bin/Rscript classic_aa.R');
% [status, msg] = system('Rscript classic_aa.R');
% if status ~= 0
%     disp(msg);
% end
!/opt/local/bin/Rscript classic_aa_test.R
load classic_aa_test_Output.mat matLatSam
!rm classic_aa_test_Input.mat classic_aa_test_Output.mat