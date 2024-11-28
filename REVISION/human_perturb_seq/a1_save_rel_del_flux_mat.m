load('output/relDeltaFluxMat.mat')

load('input/human_model/Human-GEM_cobra.mat'); % we found ihuman doesnt have csense field - so we use cobra 
tgtGenes = setdiff(model.genes,{'NA','ND','TBD','Unknown'});


% core energy obj
coreEnergyObj = {'MAR03964'};

% total lipid objectives
lipidObj = {'DMN_MAM10014c[c]'};

% total energy objectives
energyObj = coreEnergyObj;

% total protein 
proteinSynObj = {'DMN_MAM10013c[c]'};

% total ECM 
ECM_obj = readtable('input\human_model\ECM_mets.xlsx','ReadVariableNames',true);
ECMObj = ECM_obj.EX_rxns';

% total nucleic acid
nucleicAcidObj = {'DMN_MAM01721n[n]','DMN_MAM02847c[c]'};


allObjs = unique([lipidObj, energyObj, proteinSynObj,ECMObj,nucleicAcidObj]);

%% calculate the merged obj values for each obj
% merge the subobjectives (max delta flux)
rdf_coreEnergy = max(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)),[],2);
% print number of sig genes
n_coreEnergy = sum(rdf_coreEnergy > 1e-3)


% merge the subobjectives (max delta flux)
rdf_lipidObj = max(relDeltaFluxMat(:,ismember(allObjs, lipidObj)),[],2);
% print number of sig genes
n_lipidObj = sum(rdf_lipidObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_energyObj = max(relDeltaFluxMat(:,ismember(allObjs, energyObj)),[],2);
% print number of sig genes
n_energyObj = sum(rdf_energyObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_proteinModificationObj = max(relDeltaFluxMat(:,ismember(allObjs, ECMObj)),[],2);
% print number of sig genes
n_proteinModificationObj = sum(rdf_proteinModificationObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_proteinSynObj = max(relDeltaFluxMat(:,ismember(allObjs, proteinSynObj)),[],2);
% print number of sig genes
n_proteinSynObj = sum(rdf_proteinSynObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_nucleicAcidObj = max(relDeltaFluxMat(:,ismember(allObjs, nucleicAcidObj)),[],2);
% print number of sig genes
n_nucleicAcidObj = sum(rdf_nucleicAcidObj > 1e-3)

%% generate the output matrix 
t = table(rdf_energyObj, rdf_lipidObj, rdf_proteinSynObj, rdf_nucleicAcidObj, rdf_proteinModificationObj);
t = table2array(t);
t(t<1e-9) = 0;
t = array2table(t);
t.Properties.VariableNames = ["energy","lipid" ,"pro_syn", "nucl_acid","ECM"];
t.Properties.RowNames = tgtGenes;
writetable(t, 'output/delta_flux_scoreMat.csv','WriteRowNames',true)

%% save polymerases
pols = [model.genes(model.rxnGeneMat(strcmp(model.rxns,'MAR07160'),:)==1);
    model.genes(model.rxnGeneMat(strcmp(model.rxns,'MAR07161'),:)==1);
     model.genes(model.rxnGeneMat(strcmp(model.rxns,'MAR13086'),:)==1)];% RNA transport
writetable(array2table(pols),'output/DNA_RNA_pol.csv','WriteRowNames',true)





