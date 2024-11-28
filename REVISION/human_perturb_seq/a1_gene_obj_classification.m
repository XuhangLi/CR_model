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

%% notes

% the generic metabolic network model may not apply to cancer cell line. a
% better test of CR model will be using the cellline-specific model of K562
% cells and rerun the FBA. Currently the core function assignment might be
% compromised, for example, some genes may be not important or may be
% important in the actual cancer network, distinct from the generic
% network.

%% SETUP SOLVER
% human model is numerically challenging, so we set the solver stringency
% to highest; the default (1e-7) is known to show problem for optimizing
% some complex glycans  
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
%%
% Load model
load('input/human_model/Human-GEM_cobra.mat'); % we found ihuman doesnt have csense field - so we use cobra 

model = addDefaultConstraint(model,-0.001, -0.005); % arbiturary flux <1% glucose carbon

% the media constraints have a lot to optimize; likely need some real data
% guidance (like NCI60) to be coherent with the typical metabolism such as
% glycolysis tca and etc; now the huge non-glucose flux dominates the
% energy source 



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
allPrecursors = model.mets(model.S(:, strcmp(model.rxns, 'MAR13082'))<0);
allPrecursors = setdiff(allPrecursors,{'MAM01371c','MAM02040c',''}); % remove water, atp
for i = 1:length(allPrecursors)
    model = addReaction(model,['DMN_',allPrecursors{i}],'reactionFormula',[allPrecursors{i},' ->'],'geneRule', 'NA','printLevel',1);
end


% define the objective set 

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

% % total cofactors 
% cofactorObj = {'DMN_MAM10012c[c]'};
% 
% % total metabolite pool 
% metaboliteObj = {'DMN_MAM10015c[c]'};


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
parpool(12)
allObjs = unique([lipidObj, energyObj, proteinSynObj,ECMObj,nucleicAcidObj]);
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


% %% compare with reference set 
% objClass = {'lipid'};
% refType = 'coexpression_cluster';
% rel_delta_flux = rdf_lipidObj;
% 
% reflabels = readtable('input/benchmark_obj_labels.csv'); % gene name may have some mismatches - fix soon
% geneID = readtable('./../../input_data/otherTbls/iCEL_IDtbl.csv');
% [A B] = ismember(reflabels.WBID, geneID.WormBase_Gene_ID);
% reflabels.gene_name(A) = geneID.ICELgene(B(A));
% myset = reflabels.gene_name(ismember(reflabels.label,objClass) & strcmp(reflabels.reference, refType));
% 
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

allObjs = unique([lipidObj, energyObj, proteinSynObj,ECMObj,nucleicAcidObj]);

% core energy genes (only energy genes, to be excluded from other labels)
obj_energy_genes_core = unique([tgtGenes(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)) == 1);... # principle #1
                            tgtGenes(relDeltaFluxMat(:,ismember(allObjs, coreEnergyObj)) > max(relDeltaFluxMat(:,ismember(allObjs, [lipidObj, proteinSynObj,ECMObj,nucleicAcidObj])),[],2)); ... # principle #2
                            ]);

% core protein syn (only pro_syn, to be excluded from other labels)
tmp = model;
tmp.subSystems = [tmp.subSystems{:}]';
obj_proteinSyn_genes_core = model.genes(any(model.rxnGeneMat(strcmp(tmp.subSystems,'Aminoacyl-tRNA biosynthesis'),:)==1,1));


% define each obj genes
obj_energy_genes = tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, energyObj)),[],2) > cutoff);
obj_lipid_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, lipidObj)),[],2) > cutoff),[obj_energy_genes_core]);
obj_proteinModification_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, ECMObj)),[],2) > cutoff),[obj_energy_genes_core;obj_proteinSyn_genes_core]);
obj_proteinSyn_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, proteinSynObj)),[],2) > cutoff),obj_energy_genes_core);
obj_nucleicAcid_genes = setdiff(tgtGenes(max(relDeltaFluxMat(:,ismember(allObjs, nucleicAcidObj)),[],2) > cutoff),[obj_energy_genes_core]);

save('output/classification_genes.mat',"obj_energy_genes","obj_lipid_genes","obj_proteinSyn_genes","obj_nucleicAcid_genes","obj_proteinModification_genes")

% %% compare with reference set 
% objClass = {'lipid'};
% refType = 'coexpression_cluster';
% FBA_labels = obj_lipid_genes;
% 
% reflabels = readtable('input/benchmark_obj_labels.csv'); % gene name may have some mismatches - fix soon
% geneID = readtable('./../../input_data/otherTbls/iCEL_IDtbl.csv');
% [A B] = ismember(reflabels.WBID, geneID.WormBase_Gene_ID);
% reflabels.gene_name(A) = geneID.ICELgene(B(A));
% myset = reflabels.gene_name(ismember(reflabels.label,objClass) & strcmp(reflabels.reference, refType));
% 
% predicted = intersect(myset, FBA_labels)
% failed_to_predict = setdiff(myset, FBA_labels)
% %% check classification
% addpath tools/
% setListData = { toFactor(obj_energy_genes, tgtGenes),...
%                 toFactor(obj_redox_genes, tgtGenes),...
%                 toFactor(obj_lipid_genes, tgtGenes)
%                 };
% figure
% h1 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["energy","redox","lipid"]);
% h1.ShowIntersectionCounts = true;
% 
% figure
% setListData = { toFactor(obj_energy_genes, tgtGenes),...
%                 toFactor(obj_proteinModification_genes, tgtGenes),...
%                 toFactor(obj_nucleicAcid_genes, tgtGenes)
%                 };
% h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["energy", "proModification","nuclAcid"]);
% h2.ShowIntersectionCounts = true;
% 
% figure
% setListData = { toFactor(obj_energy_genes, tgtGenes),...
%                 toFactor(obj_proteinSyn_genes, tgtGenes),...
%                 toFactor(obj_proteinMito_genes, tgtGenes)
%                 };
% h3 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["energy","proSyn","proMito"]);
% h3.ShowIntersectionCounts = true;
% 
% figure
% setListData = { toFactor(obj_lipid_genes, tgtGenes),...
%                 toFactor(obj_proteinSyn_genes, tgtGenes),...
%                 toFactor(obj_proteinMito_genes, tgtGenes)
%                 };
% h4 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["lipid","proSyn","proMito"]);
% h4.ShowIntersectionCounts = true;
% 
% 
% figure
% setListData = { toFactor(obj_lipid_genes, tgtGenes),...
%                 toFactor(obj_proteinModification_genes, tgtGenes),...
%                 toFactor(obj_nucleicAcid_genes, tgtGenes)
%                 };
% h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["lipid", "proModification","nuclAcid"]);
% h2.ShowIntersectionCounts = true;
% 
% 
% 
% figure
% setListData = { toFactor(obj_proteinSyn_genes, tgtGenes),...
%                 toFactor(obj_proteinMito_genes, tgtGenes),...
%                 toFactor(obj_proteinModification_genes, tgtGenes),...
%                 toFactor(obj_nucleicAcid_genes, tgtGenes)
%                 };
% h2 = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
%     ["proSyn", "proMito","proModification","nuclAcid"]);
% h2.ShowIntersectionCounts = true;
%% classification matrix 
classMat = zeros(length(tgtGenes), 5);
classMat(:,1) = ismember(tgtGenes, obj_energy_genes);
classMat(:,2) = ismember(tgtGenes, obj_lipid_genes);
classMat(:,3) = ismember(tgtGenes, obj_proteinSyn_genes);
classMat(:,4) = ismember(tgtGenes, obj_nucleicAcid_genes);
classMat(:,5) = ismember(tgtGenes, obj_proteinModification_genes);

classMat = array2table(classMat);
classMat.Properties.VariableNames = {'energy','lipid','pro_syn','nucl_acid','ECM'};
classMat.Properties.RowNames = tgtGenes;


writetable(classMat, 'output/FBA_classification_matrix.csv','WriteRowNames',true)
