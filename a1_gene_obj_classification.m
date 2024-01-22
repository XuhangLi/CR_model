%% summary
% this simulation is primarily to classfy genes into each core function; but
% it may also help explain the regulatory hierachy by its multiple
% objectives classification;
% regardless of explaining hierachy, we classify genes as follows: 
% 1. any gene that has an rel_del_flux of 1 for core energy (RCC0005) as
% solely energy obj;
% 2. any gene that has an rel_del_flux for core energy (RCC0005) higher
% than any other non-energy obj as solely energy obj;
% 3. any gene that has an sig energy obj as energy
% 4. any gene that has an sig other obj as the other obj

%  expectation: most genes should be single obj; multiple obj should
%  include PPP genes and malvalonate genes


% to explain the hirachy of a gene, the simulation matrix needs to be
% revisited to find what obj has been influenced 

%% next step
% this simulation only simulates the energy catabolism and biomass
% anabolism; the degradatons of the objectives are not simulated; to do the
% latter, we need to reverse (or turn off?) all exchange reactions and only
% allow secrection of selected metabolites (basic units in degradation
% objective) and then calculate the maximal uptake of a given biomass obj
% and its suboptimal minimal flux and delta flux etc; 

%%
% Load model
load('input/iCEL1314.mat');
model = addDefaultConstraint(model,'minimalExchange@1');
% ATPm is an important constraint: high ATPm will impose strong energy
% coupling where reactions influcening energy production will be coupled
% stronger with the major energy comsumping reactions (by making them
% under energy limitation)
% maintance will also cause a background deltaflux influence between energy
% genes and any tested obj; we need to think about it!
% try remove the NGAM and bacterial digestion energy
%model = changeRxnBounds(model,'RCC0005',0,'l'); 
%model.S(ismember(model.mets,{'atp[c]','h2o[c]','adp[c]','h[c]','pi[c]'}), strcmp('DGR0007',model.rxns)) = 0;% remove the energy cost for bacteria digestion
% --> when these energy were ignored, met/sam degrdation of met will be
% called as energy as it supply cys-l that was used in iron-sulfur cluster
% formation; this seems a false positive; we either get rid of iron-sulfer
% or stoichemictrially couple it with energy in a single obj to weaken the
% dependecy on cys-l; or keep the bg energy demand as a context (doesnt
% seem smart tho)
% --> decide to not change; energy cost is universal and is a reasonable
% context; we should remove energy context as a post-processing

% add drains to uncouple biomass binding 
allPrecursors = model.mets(model.S(:, strcmp(model.rxns, 'BIO0107'))<0);
for i = 1:length(allPrecursors)
    model = addReaction(model,['DMN_',allPrecursors{i}],'reactionFormula',[allPrecursors{i},' ->'],'geneRule', 'NA','printLevel',1);
end

% objSet = strcat('DMN_',allPrecursors);
% objSet = [objSet; model.rxns(contains(model.rxns,'BIO')); {'RCC0005';'OBJ_nadph'}]; % will be changed to all reactions

% define the objective set 

% core energy obj
coreEnergyObj = {'RCC0005'};

% core redox obj
% add redox demand obj
model = addReaction(model,'OBJ_nadph','reactionFormula','nadph[c] -> nadp[c]','geneRule', 'NA','printLevel',1);
coreRedoxObj = {'OBJ_nadph'};

% total lipid objectives
model = addReaction(model,'DMN_pail45p','reactionFormula','pail45p[c] ->','geneRule', 'NA','printLevel',1);
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
%'DMN_pail45p',... % phosphatidylinositol --> which one to include? maybe all?
% DMN_pail45p doesnt work since alternative pathway has no flux loss
% we may skip the inositols as they are mostly signaling 
% but they seems behave like lipid genes and may account for up to 10% of
% DEG in lipid RNAi; we may include all PIPs.

% total energy objectives
model = addReaction(model,'DMN_ISclustprot','reactionFormula','ISclustprot[m] -> apoISclustprot[m]','geneRule', 'NA','printLevel',1);
model = addReaction(model,'DMN_focytC','reactionFormula','focytC[m] -> apocytc[m]','geneRule', 'NA','printLevel',1);
energyObj = {'RCC0005',... % energy demand
             'DMN_ISclustprot',... % iron-sulfur cofactors cannot be addressed in FBA without a demand
             'DMN_focytC',... % Ferrocytochrome C cofactor cannot be addressed in FBA without a demand
             'DMN0006',... % hemeA, cofator in cytc protein
            };
%'BIO0001',... % mito protein ==> the mito tRNA is to be
% manually labeled; protein syn is too complicated and
% introduce too many false positive (unexpected connections)
% B12? adocbl

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


% other notes:
% M153.1 and alh-13 were inconsistent with model as it is pro syn genes but 
% appears in enegry RNAi cluster; may be toxicity (accumulation of glu??)
% Y52B11C.1 and F54D5.12 missed for incomplete reconstructions
% alh-7 cannot be captured by current method; it is a shunt to convert akg
% to succ instead of TCA cycle; (akg --> sucsal --> succ), also recycle
% GABA; so it is not a sufficient energy and may related to GABA regulation
% idh-2 cannot be capatured with current model as it produce nadph instead
% of nadh, so the other route (idha) were chosen by the model
% ROS related genes are also missed as ROS is not in the stoichemotry
% system
% physiologically meaningful but circular flux such as pyc-1 cannot be
% captured as it is not flux efficient 

% for the lumped obj, we either loop through and collect hits for each (try first, so no interaction considered and no stoichemitry needed)
% , or we add them up, and constrain the maximum of each to 90% (then an
% interaction is enforced; ideally we should consider stoichemtry)

% load('control_obj_values.mat');
% load('deletion_obj_values.mat');
%%
% parpool(12)
allObjs = unique([coreEnergyObj,coreRedoxObj,lipidObj,energyObj,proteinModificationObj,proteinSynObj,proteinMitoObj,nucleicAcidObj]);
%% simulate the most efficient flux for each obj and calculate the delta flux sensitivity 
myObjs = allObjs;
suboptimality_cutoff = 0.9; % a little flux (simulating delta flux) or 90% simulating suboptimal
GPRmethod = 'association';
tgtGenes = setdiff(model.genes,{'NA','ND','TBD','Unknown'});

relDeltaFluxMat = nan(length(tgtGenes), length(myObjs));
for zz = 1:length(myObjs)
    myObj = myObjs{zz};

    % calculate the most efficient flux 
    test = changeObjective(model,myObj);
    sol = optimizeCbModel(test);
    optimal_obj = sol.f * suboptimality_cutoff; % we use the principle of suboptimality 
    test.lb(strcmp(test.rxns,myObj)) = optimal_obj;
    opt_flux = minimizeModelFlux_XL(test);
    optimal_totalFlux = sum(abs(opt_flux));
    
    environment = getEnvironment();
    % simulate the delta flux 
    rel_delta_flux = [];
    fprintf(['\n' repmat('.',1,length(tgtGenes)) '\n\n']);
    parfor i = 1:length(tgtGenes)
        restoreEnvironment(environment);
        tgtGene = tgtGenes(i);
        test = changeObjective(model,myObj);
        test.lb(strcmp(test.rxns,myObj)) = optimal_obj;
        if strcmp(GPRmethod,'regular')
            test = deleteModelGenes(test,tgtGene);
        elseif strcmp(GPRmethod,'association')
            assoRxns = model.rxns(model.rxnGeneMat(:,strcmp(model.genes,tgtGene)) == 1);
            test = changeRxnBounds(test,assoRxns,0,'b');
        end
    
        flux = minimizeModelFlux_XL(test);
        if isstruct(flux)
            rel_delta_flux(i) = 1;
        else
            del_totalFlux = sum(abs(flux));
            rel_delta_flux(i) = (del_totalFlux - optimal_totalFlux) / optimal_totalFlux;
        end
        fprintf('\b|\n');%for simple progress monitor
    end

    relDeltaFluxMat(:,zz) = rel_delta_flux;
end
save('output/relDeltaFluxMat.mat',"relDeltaFluxMat")
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

%% compare with reference set 
objClass = {'lipid'};
refType = 'coexpression_cluster';
rel_delta_flux = rdf_lipidObj;

reflabels = readtable('input/benchmark_obj_labels.csv'); % gene name may have some mismatches - fix soon
geneID = readtable('./../../input_data/otherTbls/iCEL_IDtbl.csv');
[A B] = ismember(reflabels.WBID, geneID.WormBase_Gene_ID);
reflabels.gene_name(A) = geneID.ICELgene(B(A));
myset = reflabels.gene_name(ismember(reflabels.label,objClass) & strcmp(reflabels.reference, refType));


figure
hold on
histogram(log10(1e-8+rel_delta_flux(ismember(tgtGenes, myset))),'Normalization','probability',NumBins=30)
histogram(log10(1e-8+rel_delta_flux),'Normalization','probability',NumBins=30)
hold off
legend({'subset','all'})

% based on lipid and energy, the delta flux method looks better with a
% common cutoff of 0.1% (0.001) for significant compromise of flux
% efficiency; the other parameter, suboptimality cutoff is generally very
% robust, eg 0.5-0.9;

%% calculate the final obj classifications
% 1. any gene that has an rel_del_flux of 1 for core energy (RCC0005) as
% solely energy obj;
% 2. any gene that has an rel_del_flux for core energy (RCC0005) higher
% than any other non-energy obj as solely energy obj;
% 3. any gene that has an sig energy obj as energy
% 4. any gene that has an sig other obj as the other obj
% 5. all tRNA synthase only assigned to pro_syn (to avoid their appearence
% in all protein-dependent processes like collagen)

%  expectation: most genes should be single obj; multiple obj should
%  include PPP genes and malvalonate genes
cutoff = 1e-3;

allObjs = unique([coreEnergyObj,coreRedoxObj,lipidObj,energyObj,proteinModificationObj,proteinSynObj,proteinMitoObj,nucleicAcidObj]);

% core energy genes (only energy genes, to be excluded from other labels)
obj_energy_genes_core = unique([tgtGenes(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)) == 1);... # principle #1
                            tgtGenes(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)) > max(relDeltaFluxMat(:,ismember(allObjs, [lipidObj,proteinModificationObj,proteinSynObj,proteinMitoObj,nucleicAcidObj])),[],2)); ... # principle #2
                            ]);
% core protein syn (only pro_syn, to be excluded from other labels)
tmp = model;
tmp.subSystems = [tmp.subSystems{:}]';
obj_proteinSyn_genes_core = model.genes(any(model.rxnGeneMat(strcmp(tmp.subSystems,'Aminoacyl-tRNA biosynthesis'),:)==1,1));


% define each obj genes
obj_energy_genes = tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, energyObj)),[],2) > cutoff);
obj_redox_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, coreRedoxObj)),[],2) > cutoff),[obj_energy_genes_core;obj_proteinSyn_genes_core]);
obj_lipid_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, lipidObj)),[],2) > cutoff),[obj_energy_genes_core;obj_proteinSyn_genes_core]);
obj_proteinModification_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, proteinModificationObj)),[],2) > cutoff),[obj_energy_genes_core;obj_proteinSyn_genes_core]);
obj_proteinSyn_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, proteinSynObj)),[],2) > cutoff),obj_energy_genes_core);
obj_proteinMito_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, proteinMitoObj)),[],2) > cutoff),obj_energy_genes_core);
obj_nucleicAcid_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, nucleicAcidObj)),[],2) > cutoff),[obj_energy_genes_core;obj_proteinSyn_genes_core]);

save('output/classification_genes.mat',"obj_energy_genes","obj_redox_genes","obj_lipid_genes","obj_proteinModification_genes","obj_proteinSyn_genes","obj_proteinMito_genes","obj_nucleicAcid_genes")

%% compare with reference set 
objClass = {'lipid'};
refType = 'coexpression_cluster';
FBA_labels = obj_lipid_genes;

reflabels = readtable('input/benchmark_obj_labels.csv'); % gene name may have some mismatches - fix soon
geneID = readtable('./../../input_data/otherTbls/iCEL_IDtbl.csv');
[A B] = ismember(reflabels.WBID, geneID.WormBase_Gene_ID);
reflabels.gene_name(A) = geneID.ICELgene(B(A));
myset = reflabels.gene_name(ismember(reflabels.label,objClass) & strcmp(reflabels.reference, refType));

predicted = intersect(myset, FBA_labels)
failed_to_predict = setdiff(myset, FBA_labels)
%% check classification
addpath tools/
setListData = { toFactor(obj_energy_genes, tgtGenes),...
                toFactor(obj_redox_genes, tgtGenes),...
                toFactor(obj_lipid_genes, tgtGenes)
                };
figure
h1 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["energy","redox","lipid"]);
h1.ShowIntersectionCounts = true;

figure
setListData = { toFactor(obj_energy_genes, tgtGenes),...
                toFactor(obj_proteinModification_genes, tgtGenes),...
                toFactor(obj_nucleicAcid_genes, tgtGenes)
                };
h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["energy", "proModification","nuclAcid"]);
h2.ShowIntersectionCounts = true;

figure
setListData = { toFactor(obj_energy_genes, tgtGenes),...
                toFactor(obj_proteinSyn_genes, tgtGenes),...
                toFactor(obj_proteinMito_genes, tgtGenes)
                };
h3 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["energy","proSyn","proMito"]);
h3.ShowIntersectionCounts = true;

figure
setListData = { toFactor(obj_lipid_genes, tgtGenes),...
                toFactor(obj_proteinSyn_genes, tgtGenes),...
                toFactor(obj_proteinMito_genes, tgtGenes)
                };
h4 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["lipid","proSyn","proMito"]);
h4.ShowIntersectionCounts = true;


figure
setListData = { toFactor(obj_lipid_genes, tgtGenes),...
                toFactor(obj_proteinModification_genes, tgtGenes),...
                toFactor(obj_nucleicAcid_genes, tgtGenes)
                };
h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["lipid", "proModification","nuclAcid"]);
h2.ShowIntersectionCounts = true;



figure
setListData = { toFactor(obj_proteinSyn_genes, tgtGenes),...
                toFactor(obj_proteinMito_genes, tgtGenes),...
                toFactor(obj_proteinModification_genes, tgtGenes),...
                toFactor(obj_nucleicAcid_genes, tgtGenes)
                };
h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    ["proSyn", "proMito","proModification","nuclAcid"]);
h2.ShowIntersectionCounts = true;
%% classification matrix 
classMat = zeros(length(tgtGenes), 7);
classMat(:,1) = ismember(tgtGenes, obj_energy_genes);
classMat(:,2) = ismember(tgtGenes, obj_redox_genes);
classMat(:,3) = ismember(tgtGenes, obj_lipid_genes);
classMat(:,4) = ismember(tgtGenes, obj_proteinModification_genes);
classMat(:,5) = ismember(tgtGenes, obj_proteinSyn_genes);
classMat(:,6) = ismember(tgtGenes, obj_proteinMito_genes);
classMat(:,7) = ismember(tgtGenes, obj_nucleicAcid_genes);

classMat = array2table(classMat);
classMat.Properties.VariableNames = {'energy','nadph','lipid','pro_modi','pro_syn','pro_mito','nucl_acid'};
classMat.Properties.RowNames = tgtGenes;

writetable(classMat, 'output/FBA_classification_matrix.csv','WriteRowNames',true)
%% inspection
setdiff(myset, tgtGenes(rel_delta_flux > 1e-3))
%%
intersect(myset, tgtGenes(rel_delta_flux > 1e-3))

%% note
% the lipid obj simulation will include genes for lipid (incl. redox and 
% cofactors) and energy; but to fit the model, lipid should be defined as
% only inflencing lipid, therefore, overlap genes with energy should be
% substract. In this case, the genes in energy simulation shouldnt be in
% the lipid cluster significantly (esp. in RNAi_cluster lipid) ==> it is
% unbelievably true!!!(now for RNAi_cluster lipid)(and as expected, PPP
% genes are exceptions (truely both in lipid and energy); hacd-1 and fah-1
% are also reasonable exceptions;); ==> the obj classification of exception
% genes can be mannually adjusted at the end

% from FBA we realized that the up-regulation of AA deg, especially serine
% and glycine deg, may be for compensating nadph; as they are major ways to
% support nadph production by folate cycle (either through making cytosolic
% mlthf directly or making mitochondrial nadph). (serine syn, folate nadph
% syn, and some case for TCA nadph syn (idh-1), are indeed up in lipid
% RNAi)

% from FBA we realized that a large amout of nadph could be from acald
% oxidation that can be produced from aprop (ala+nh4); this reaction was 
% catalyzed by nit-1, whcih was up sig in many lipid RNA that we dont
% underatand before


% quite a lot genes involved in NTP-NDP/NMP conversions were found in
% energy RNAi cluster but no idea why connected to energy as well as not
% captured in FBA

%%
listRxn(test,opt_flux,'nadph[c]')
%%
%% old codes

%%
optimizeCbModel(test)
optimizeCbModel(test)
test = deleteModelGenes(test, 'RM00248')
test = deleteModelGenes(test, 'gdh-1')
optimizeCbModel(test)
sol = optimizeCbModel(test)
listRxn(test,sol.full,'nadph[c]')
listRxn(test,sol.full,'4hphac[c]')
listRxn(test,sol.full,'nadph[c]')
%%
listRxn(test,opt_flux,'nadph[c]')

% test = deleteModelGenes(test, 'K07E3.4') % folate 1
% test = deleteModelGenes(test, 'alh-3') % folate 2
% test = deleteModelGenes(test, 'aco-2')
%test = deleteModelGenes(test, 'gldc-1')
%test = changeRxnBounds(test,'BIO0010',0.7,'l'); % default 10x bacterial uptake causes ETC infeasible, we lower down to 5x
%test = changeRxnBounds(test,'RM00751',0,'b'); % default 10x bacterial uptake causes ETC infeasible, we lower down to 5x
%test = changeRxnBounds(test,{'RM03293','RM00248','RC00711','RM00711','RM00216'},0,'b'); % default 10x bacterial uptake causes ETC infeasible, we lower down to 5x

%% 
relVm = Vmax_mat ./ repmat(Vmax_control, size(Vmax_mat,1),1);

Vmax_check = array2table(relVm);
Vmax_check.Properties.VariableNames = objSet;
Vmax_check.Properties.RowNames = geneSet;
%% code backup for single obj analysis
% %% simulate the most efficient flux for each obj
% myObj = 'RCC0005';
% suboptimality_cutoff = 0.9; % a little flux (simulating delta flux) or 90% simulating suboptimal
% sigFluxThreshold = 0.1; % 10% of system influx (bacteria) --> this is a tunable parameter 
% 
% test = changeObjective(model,myObj);
% % test = deleteModelGenes(test, 'bckd-1A');
% sol = optimizeCbModel(test);
% optimal_obj = sol.obj * suboptimality_cutoff; % we use the principle of suboptimality 
% test.lb(strcmp(test.rxns,myObj)) = optimal_obj;
% opt_flux = minimizeModelFlux_XL(test);
% 
% optimal_totalFlux = sum(abs(opt_flux))
% % sum(any(model.rxnGeneMat(abs(opt_flux)> abs(opt_flux(strcmp(model.rxns,'EXC0050'))) * sigFluxThreshold,:)==1,1))
% 
% %% calculate the total flux for each gene 
% gene_flux = [];
% for i = 1:length(model.genes)
%     myrxns = model.rxns(model.rxnGeneMat(:,i)==1);
%     gene_flux(i,1) = sum(abs(opt_flux(ismember(model.rxns,myrxns))));
% end
% 
% n_total = sum(gene_flux > abs(opt_flux(strcmp(model.rxns,'EXC0050'))) * sigFluxThreshold)
% %% compare with reference set 
% objClass = {'energy'};
% refType = 'RNAi_cluster';
% 
% reflabels = readtable('benchmark_obj_labels.csv'); % gene name may have some mismatches - fix soon
% geneID = readtable('./../../../../input_data/otherTbls/iCEL_IDtbl.csv');
% [A B] = ismember(reflabels.WBID, geneID.WormBase_Gene_ID);
% reflabels.gene_name(A) = geneID.ICELgene(B(A));
% myset = reflabels.gene_name(ismember(reflabels.label,objClass) & strcmp(reflabels.reference, refType));
% 
% figure
% hold on
% histogram(log10(1e-8+gene_flux(ismember(model.genes, myset))),'Normalization','probability',NumBins=30)
% histogram(log10(1e-8+gene_flux),'Normalization','probability',NumBins=30)
% hold off
% legend({'subset','all'})
% %% test the delta flux 
% % simple delta or MOMA delta?? (simple delta may have the infeasible issue
% % but MOMA delta will always have a number) ==> MOMA seems not very solid
% % as we are not interested in how much the overfitted minimal flux
% % distribution is disrupted, but how much cost of flux has to be paid to
% % rewire 
% 
% tgtGenes = setdiff(model.genes,{'NA','ND','TBD','Unknown'});
% GPRmethod = 'association';
% 
% rel_delta_flux = [];
% for i = 1:length(tgtGenes)
%     tgtGene = tgtGenes(i);
%     test = changeObjective(model,myObj);
%     test.lb(strcmp(test.rxns,myObj)) = optimal_obj;
%     if strcmp(GPRmethod,'regular')
%         test = deleteModelGenes(test,tgtGene);
%     elseif strcmp(GPRmethod,'association')
%         assoRxns = model.rxns(model.rxnGeneMat(:,strcmp(model.genes,tgtGene)) == 1);
%         test = changeRxnBounds(test,assoRxns,0,'b');
%     end
% 
%     flux = minimizeModelFlux_XL(test);
%     if isstruct(flux)
%         rel_delta_flux(i) = 1;
%     else
%         del_totalFlux = sum(abs(flux));
%         rel_delta_flux(i) = (del_totalFlux - optimal_totalFlux) / optimal_totalFlux;
%     end
%     i/length(tgtGenes)
% end
% 
% sum(rel_delta_flux > 1e-3)
% %% compare with reference set 
% 
% figure
% hold on
% histogram(log10(1e-8+rel_delta_flux(ismember(tgtGenes, myset))),'Normalization','probability',NumBins=30)
% histogram(log10(1e-8+rel_delta_flux),'Normalization','probability',NumBins=30)
% hold off
% legend({'subset','all'})
% 
% % based on lipid and energy, the delta flux method looks better with a
% % common cutoff of 0.1% (0.001) for significant compromise of flux
% % efficiency; the other parameter, suboptimality cutoff is generally very
% % robust, eg 0.5-0.9;
