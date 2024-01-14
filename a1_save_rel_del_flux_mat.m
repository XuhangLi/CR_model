load('output/relDeltaFluxMat.mat')

load('input/iCEL1314.mat');
tgtGenes = setdiff(model.genes,{'NA','ND','TBD','Unknown'});

% core energy obj
coreEnergyObj = {'RCC0005'};

% core redox obj
% add redox demand obj
coreRedoxObj = {'OBJ_nadph'};

% total lipid objectives
lipidObj = {'DMN_PhosphoL[c]',...% phospholipids
            'EXC0213','EXC0207','EXC0211'... % sphingolipids --> include several major ones (glucosphinglipid, sphingosine, sphingoline, shingomynine already in biomass)
                               ... % alh-4 and spl-1 are around sph1p_ce,
                               ... % however, the reconstruction look weird so
                               ... % there is no simple fix to rescue these
                               ... % two in the RNAi_cluster; 
            'EX04317' ... % ether lipid second branch, the product of the other branch (alkenac2gpe) has been in phospholipid biomass precursor
            'DMN_tag[c]' ... % making TAG is also considered as lipid making
            'DMN0044',... % cyclic cofactors (COA) cannot be addressed in FBA without a demand
            };

% total energy objectives
energyObj = {'RCC0005',... % energy demand
             'DMN_ISclustprot',... % iron-sulfur cofactors cannot be addressed in FBA without a demand
             'DMN_focytC',... % Ferrocytochrome C cofactor cannot be addressed in FBA without a demand
             'DMN0006',... % hemeA, cofator in cytc protein
            };

% total modified protein
proteinModificationObj = {  'DMN_Glycans[c]',... % all glycans
                            'DMN0075', ... % cofactor dolp cannot be addressed without a demand
                                       ... % we couldnt address the cyclic intermidiate, UDP-sugar/GDP-sugar etc; as a result, hprt-1 was missed here.
                                       ... % we avoid overfitting it for now, since we didnt systematically address all cyclic intermediates other than cofactors
                                       ... % we may choose to address udpg if it is linked to too many false negatives
                            'DMN_collg[c]'}; % maybe collg should be out

% total protein 
proteinSynObj = {'DMN_prot_other[c]'};

% total mito protein
proteinMitoObj = {'DMN_prot_mito[m]'};

% total nucleic acid
nucleicAcidObj = {'DMN_RNA[c]','DMN_DNA[c]'};


allObjs = unique([coreEnergyObj,coreRedoxObj,lipidObj,energyObj,proteinModificationObj,proteinSynObj,proteinMitoObj,nucleicAcidObj]);

%% calculate the merged obj values for each obj
% merge the subobjectives (max delta flux)
rdf_coreEnergy = max(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)),[],2);
% print number of sig genes
n_coreEnergy = sum(rdf_coreEnergy > 1e-3)

% merge the subobjectives (max delta flux)
rdf_coreRedox = max(relDeltaFluxMat(:,ismember(allObjs, coreRedoxObj)),[],2);
% print number of sig genes
n_coreRedox = sum(rdf_coreRedox > 1e-3)

% merge the subobjectives (max delta flux)
rdf_lipidObj = max(relDeltaFluxMat(:,ismember(allObjs, lipidObj)),[],2);
% print number of sig genes
n_lipidObj = sum(rdf_lipidObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_energyObj = max(relDeltaFluxMat(:,ismember(allObjs, energyObj)),[],2);
% print number of sig genes
n_energyObj = sum(rdf_energyObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_proteinModificationObj = max(relDeltaFluxMat(:,ismember(allObjs, proteinModificationObj)),[],2);
% print number of sig genes
n_proteinModificationObj = sum(rdf_proteinModificationObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_proteinSynObj = max(relDeltaFluxMat(:,ismember(allObjs, proteinSynObj)),[],2);
% print number of sig genes
n_proteinSynObj = sum(rdf_proteinSynObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_proteinMitoObj = max(relDeltaFluxMat(:,ismember(allObjs, proteinMitoObj)),[],2);
% print number of sig genes
n_proteinMitoObj = sum(rdf_proteinMitoObj > 1e-3)

% merge the subobjectives (max delta flux)
rdf_nucleicAcidObj = max(relDeltaFluxMat(:,ismember(allObjs, nucleicAcidObj)),[],2);
% print number of sig genes
n_nucleicAcidObj = sum(rdf_nucleicAcidObj > 1e-3)

%% generate the output matrix 
t = table(rdf_energyObj, rdf_lipidObj, rdf_proteinModificationObj, rdf_proteinSynObj, rdf_nucleicAcidObj);
t = table2array(t);
t(t<1e-9) = 0;
t = array2table(t);
t.Properties.VariableNames = ["energy","lipid" ,"pro_modi","pro_syn", "nucl_acid"];
t.Properties.RowNames = tgtGenes;
writetable(t, 'output/delta_flux_scoreMat.csv','WriteRowNames',true)



