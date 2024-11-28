
load('Human-GEM.mat');
model = ravenCobraWrapper(ihuman);
save("Human-GEM_cobra.mat",'model');
% note: we found that the Boundary metabolites are now removed from riven
% models and the field names are well aligned with cobra model. so we just
% go with the original model until some error occurs (becasue we dont know
% if the ravencobrawrapper function is still doing the right thing now). 