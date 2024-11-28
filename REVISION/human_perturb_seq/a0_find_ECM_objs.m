%% find the metabolites related to ECM core function 

% ECM are generally materials synthesized in ER and/or golgi and ended up
% secrected. Therefore, we can check the metabolites in common appearence
% in ER+EXT and golgi+EXT


load('input/human_model/Human-GEM_cobra.mat'); % we found ihuman doesnt have csense field - so we use cobra 

%% find metabolites

ER_met = regexprep(model.metNames(contains(model.metNames, 'Endoplasmic reticulum')),' \[Endoplasmic reticulum\]','');
Golgi_met = regexprep(model.metNames(contains(model.metNames, 'Golgi apparatus')),' \[Golgi apparatus\]','');
Ex_met = regexprep(model.metNames(contains(model.metNames, 'Extracellular')),' \[Extracellular\]','');

% candidates
candi = union(intersect(ER_met, Ex_met),intersect(Golgi_met, Ex_met));

% the ECM objective functions were manually curated to 9 representative
% metabolites, representing all secrected glycans or membrance surface. 