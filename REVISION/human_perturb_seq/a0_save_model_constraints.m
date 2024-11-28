%% save the fully constrianed model that is ready for CR simulation, with propoer notes attached 
initCobraToolbox(false);

VitLim = -0.001;
essenlipidLim = -0.005;

%% Load model
load('input/human_model/Human-GEM_cobra.mat'); % we found ihuman doesnt have csense field - so we use cobra 
model.rxnNotes = repmat({'Original human 1 reaction and default constraints.'},length(model.rxns),1);

% model = addDefaultConstraint(model,-0.001, -0.005); % arbiturary flux <1% glucose carbon
% copied codes here
model.lb(findExcRxns(model)) = 0;
model.rxnNotes(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = repmat({'Only allow secretion on exchange reactions to simulate parsimonious nutrient condition.'},sum(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)),1);


 % impose real-data flux constraints 
% NCI-60 exchange rates has the same K562 cells under a near-identical
% culture condition (complete medium containing RPMI-1640 (GIBCO) with 
% 2 mM L-glutamine (~0.3 g/L) and 5% fetal bovine serum.

% load exchange flux data from NCI-60
% in mmole/gDW/h, glucose is ~0.5-1

% load target 
load('./input/human_cell_flux_data/normalizedFlux.mat');
rxnLabel = regexprep(rxnLabel,'\(|\)|\[|\]|-','_');
IDtbl = readtable('input/human_cell_flux_data/recon2iHumanEx.xlsx');
[A B] = ismember(rxnLabel, IDtbl.Recon2ID);
rxnLabel = IDtbl.Human1_ID(B(A));

for i = 1:size(rxnLabel,1)
    myrxn = rxnLabel{i};
    % we only constrain the uptake flux for now. constrain both uptake
    % and secrection causes infeasible model, whcih might be due to
    % some blocked EX and we dont want to fix 
    if flux_lb(i, strcmp(cellLabel, 'LINE_K562')) < 0
        if any(strcmp(model.rxns,myrxn))
            model.lb(strcmp(model.rxns,myrxn)) = flux_lb(i, strcmp(cellLabel, 'LINE_K562'));
            model.rxnNotes(strcmp(model.rxns,myrxn)) = {'Experimental exchange flux (uptake rates) from PMID: 28120890.'};

        else % add sink
            model = addSinkReactions(model,[myrxn,'[c]'],flux_lb(i, strcmp(cellLabel, 'LINE_K562')),1000);
            model.rxnNotes(strcmp(model.rxns,myrxn)) = {'New sink reaction reflecting experimentally observed uptake flux in PMID: 28120890.'};

        end
    end
end

model = changeRxnBounds(model,'MAR09145',VitLim,'l'); % pantothenate   
model.rxnNotes(strcmp(model.rxns,'MAR09145')) = {'Special constraint to avoid biomass limiting by cofactor uptake flux. Although the uptake flux of this cofactor is determined in PMID: 28120890, we found that the experimental flux constraint strongly limits the biomass production, likely caused by inaccurate biomass composition of this cofactor in human1. To avoid unrealistic constraints, we used an arbiturary small flux as the constraint instead of the experimental value.'}

model = changeRxnBounds(model,'MAR09361',VitLim,'l'); % inositol 
model.rxnNotes(strcmp(model.rxns,'MAR09361')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09146',VitLim,'l'); % folate 
model.rxnNotes(strcmp(model.rxns,'MAR09146')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09144',VitLim,'l'); % pyridoxine  
model.rxnNotes(strcmp(model.rxns,'MAR09144')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09159',VitLim,'l'); % thiamine  
model.rxnNotes(strcmp(model.rxns,'MAR09159')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09109',VitLim,'l'); % biotin   
model.rxnNotes(strcmp(model.rxns,'MAR09109')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09143',VitLim,'l'); % riboflavin   
model.rxnNotes(strcmp(model.rxns,'MAR09143')) = {'Allow arbiturary small flux for essential vitamin uptake.'}
model = changeRxnBounds(model,'MAR09269',VitLim,'l'); % B12   
model.rxnNotes(strcmp(model.rxns,'MAR09269')) = {'Allow arbiturary small flux for essential vitamin uptake.'}

% all other cofactors in cofactor pool avaiable in very trace amout 
cofactors = model.metNames(model.S(:,strcmp(model.rxns,'MAR10065')) < 0);
cofactorEx = regexprep(cofactors, ' \[.+\]$',' [Extracellular]');
cofactorExRxns = intersect(model.rxns(any(model.S(ismember(model.metNames,cofactorEx),:) < 0)), model.rxns(findExcRxns(model)));
% skip ones with existing bonds (from data)
cofactorExRxns = setdiff(cofactorExRxns, model.rxns(model.lb ~= 0));
model = changeRxnBounds(model,cofactorExRxns,VitLim,'l'); 
model.rxnNotes(ismember(model.rxns,cofactorExRxns)) = {'Allow arbiturary small flux for essential vitamin/cofactor uptake.'}



% vitamins that is transformed (so not in cofactorEx) but needed 
model = changeRxnBounds(model,{'MAR09151','MAR09154'},VitLim,'l'); % vitE
model = changeRxnBounds(model,'MAR10492',VitLim,'l'); % retinol   
model.rxnNotes(ismember(model.rxns,{'MAR09151','MAR09154','MAR10492'})) = {'Allow arbiturary small flux for essential vitamin/cofactor uptake.'}


% second, allow essential FA
% they are C18 (three times carbon for glucose. To avoid significant carbon
% flux (<1% glucose) ~1.5/3/100 = 0.005;
model = changeRxnBounds(model,'MAR09035',essenlipidLim,'l'); % linoneate
model = changeRxnBounds(model,'MAR09036',essenlipidLim,'l'); % linolenate
model.rxnNotes(ismember(model.rxns,{'MAR09035','MAR09036'})) = {'Allow arbiturary small flux for essential lipids uptake.'}


% allow all inorganic nutrient to be unlimited 
model = changeRxnBounds(model,'MAR09150',-1000,'l'); % cl 
model = changeRxnBounds(model,'MAR09077',-1000,'l'); % na 
model = changeRxnBounds(model,'MAR09078',-1000,'l'); % hco3 
model = changeRxnBounds(model,'MAR09072',-1000,'l'); % pi 
model = changeRxnBounds(model,'MAR09081',-1000,'l'); % k+ 
model = changeRxnBounds(model,'MAR13072',-1000,'l'); % mg2+ 
model = changeRxnBounds(model,'MAR09074',-1000,'l'); % so42- 
model = changeRxnBounds(model,'MAR09082',-1000,'l'); % ca2+ 
model = changeRxnBounds(model,'MAR09149',-1000,'l'); % no2- 
model = changeRxnBounds(model,'MAR09047',-1000,'l'); % h2o
model = changeRxnBounds(model,'MAR09048',-1000,'l'); % o2
model = changeRxnBounds(model,'MAR09079',-1000,'l'); % h+
model.rxnNotes(ismember(model.rxns,{'MAR09150','MAR09077','MAR09078','MAR09072', ...
                                    'MAR09081','MAR13072','MAR09074','MAR09082', ...
                                    'MAR09149','MAR09047','MAR09048','MAR09079'})) = {'Allow unlimited uptake flux for inorganic ions.'}


% next a few aa are not measured, use a number from similar AA
model = changeRxnBounds(model,'MAR09038',-0.01,'l'); % histidine; % a conservative number based on the flux of other essential AA (~0.01-0.03)
model.rxnNotes(ismember(model.rxns,{'MAR09038'})) = {'Allow arbiturary moderate flux for histidine uptake, which is not measured but exist in the media. This flux is within the range of general amino acid uptake flux in the dataset (~0.01-0.03).'}
model = changeRxnBounds(model,'MAR09363',-0.01,'l'); % cystine 
model.rxnNotes(ismember(model.rxns,{'MAR09363'})) = {'Allow arbiturary moderate flux for cystine uptake, which is not measured but exist in the media. This flux is within the range of general amino acid uptake flux in the dataset (~0.01-0.03).'}
model = changeRxnBounds(model,'MAR09351',-0.01,'l'); % glutathione 
model.rxnNotes(ismember(model.rxns,{'MAR09351'})) = {'Allow arbiturary moderate flux for GSH uptake, which is not measured but exist in the media. This flux is within the range of general amino acid uptake flux in the dataset (~0.01-0.03).'}

% set obj
model = changeObjective(model,'MAR13082');% we stay with the default biomass instead of total biomass (MAR10024) for now


% add drains to uncouple biomass binding 
allPrecursors = model.mets(model.S(:, strcmp(model.rxns, 'MAR13082'))<0);
allPrecursors = setdiff(allPrecursors,{'MAM01371c','MAM02040c',''}); % remove water, atp
for i = 1:length(allPrecursors)
    model = addReaction(model,['DMN_',allPrecursors{i}],'reactionFormula',[allPrecursors{i},' ->'],'geneRule', 'NA','printLevel',1);
    model.rxnNotes(end) = {'New reaction added to uncouple the draining fluxes of different biomass precursors.'};
end


%% save
model_out = struct();
model_out.rxnID  = model.rxns;
model_out.rxnFormula  = printRxnFormula_XL(model,model.rxns,false);
model_out.LB  = model.lb;
model_out.UB  = model.ub;
model_out.note  = model.rxnNotes;
model_out = struct2table(model_out);

writetable(model_out,'output/human1_constriants.csv');
