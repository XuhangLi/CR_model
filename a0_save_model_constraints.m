%% save the fully constrianed model that is ready for CR simulation, with propoer notes attached 
initCobraToolbox(false);

%% Load model
load('input/iCEL1314.mat');
model.rxnNotes = repmat({'Original iCEL1314 reaction and default constraints.'},length(model.rxns),1);

% model = addDefaultConstraint(model,'minimalExchange@1');
% copied codes here
model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = 0;
model.rxnNotes(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = repmat({'Only allow secretion on exchange reactions to simulate parsimonious nutrient condition.'},sum(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)),1);

% model.lb(cellfun(@(x) strcmp(x(1:3),'SNK'),model.rxns)) = 0;
% sink reactions were turned off by default

model = changeObjective(model,'BIO0010'); % no need to annotation - not a constraint 

model = changeRxnBounds(model,'DMN0033',0,'l');
model = changeRxnBounds(model,'DMN0033',0.01,'u');
model.rxnNotes(strcmp(model.rxns,'DMN0033')) = {'Only allow minimal cytosolic q regeneration since this reaction is not masss-balanced and causes redox leak. It represents oxidative stress per model reconstruction note #31.'}

model = changeRxnBounds(model,'RCC0005',10,'l'); % assume NGAM as 10x bacterial uptake
model = changeRxnBounds(model,'RCC0005',1000,'u');
model.rxnNotes(strcmp(model.rxns,'RCC0005')) = {'Assume NGAM is 10x the rate of bacterial uptake (Yilmaz et al., Cell System, 2016)'}

model = changeRxnBounds(model,'RM00112',0,'l');
model = changeRxnBounds(model,'RM00112',0.01,'u'); % nnt-1 reaction; assume producing nadph and assume not massive production (since H+ gradient coupling was not reconstructed)
model.rxnNotes(strcmp(model.rxns,'RM00112')) = {'Constrain the flux through nnt-1 reaction to maintain proper redox balance. Since H+ gradient coupling was not reconstructed in the model, this reaction can convert nadh to nadph without any cost. Therefore, we constain it to minimal level of nadph production to avoid massive thermodynamically infeasible nadph production.'}

% model = changeRxnBounds(model,'RM04432',-1000,'l');
% model = changeRxnBounds(model,'RM04432',1000,'u');
% nothing is changed 

model = changeRxnBounds(model,'RMC0005',0,'u');
model = changeRxnBounds(model,'RMC0005',0,'l');
model.rxnNotes(strcmp(model.rxns,'RMC0005')) = {'Block one of the two alternative ATPase reactions, to simplify the solution space.'}

% Essential exchanges:
model = changeRxnBounds(model,'EXC0050',-1,'l');
model.rxnNotes(strcmp(model.rxns,'EXC0050')) = {'Simulate one unit of bacterial diet influx.'};
model = changeRxnBounds(model,'EX00001',-1000,'l');
model = changeRxnBounds(model,'EX00007',-1000,'l');
model = changeRxnBounds(model,'EX00009',-1000,'l');
model = changeRxnBounds(model,'EX00080',-1000,'l');
model.rxnNotes(ismember(model.rxns,{'EX00001','EX00007','EX00009','EX00080'})) = {'Allow unlimited uptake of inorganic nutrients, including h2o, o2, pi and h.'};

% allowed storage (to simulate genes related to storage ultilization)
model = changeRxnBounds(model,'SNK0012',-0.01,'l');
model = changeRxnBounds(model,'SNK0013',-0.01,'l');
model = changeRxnBounds(model,'SNK0014',-0.01,'l'); % start with 0.01 TAG and this can be lower if it causes irrealisitic carbon influx
model.rxnNotes(ismember(model.rxns,{'SNK0012','SNK0013','SNK0014'})) = {'Allow limited use of storage molecules, including trehalose, TAG and glycogen.'};


% Additional exchanges:
model = changeRxnBounds(model,'EX00089',-0.01,'l');
model = changeRxnBounds(model,'EX00121',-0.01,'l');
model = changeRxnBounds(model,'EX00132',-0.01,'l');
model = changeRxnBounds(model,'EX00133',-0.01,'l');
model = changeRxnBounds(model,'EX00185',-0.01,'l');
model = changeRxnBounds(model,'EX00187',-0.01,'l');
model = changeRxnBounds(model,'EX00208',-0.01,'l');
model = changeRxnBounds(model,'EX00243',-0.01,'l');
model = changeRxnBounds(model,'EX00257',-0.01,'l');
model = changeRxnBounds(model,'EX00402',-0.01,'l');
model = changeRxnBounds(model,'EX00486',-0.01,'l');
model = changeRxnBounds(model,'EX00523',-0.01,'l');
model = changeRxnBounds(model,'EX00535',-0.01,'l');
model = changeRxnBounds(model,'EX00708',-0.01,'l');
model = changeRxnBounds(model,'EX00777',-0.01,'l');
model = changeRxnBounds(model,'EX01613',-0.01,'l');
model = changeRxnBounds(model,'EX01801',-0.01,'l');
model = changeRxnBounds(model,'EX02094',-0.01,'l');
model = changeRxnBounds(model,'EX02415',-0.01,'l');
model = changeRxnBounds(model,'EX05127',-0.01,'l');
model = changeRxnBounds(model,'EX05402',-0.01,'l');
model = changeRxnBounds(model,'EX05422',-0.01,'l');
model = changeRxnBounds(model,'EX05697',-0.01,'l');
model = changeRxnBounds(model,'EXC0152',-0.01,'l');
model = changeRxnBounds(model,'EX01132',-0.01,'l');
model = changeRxnBounds(model,'EX18125',-0.001,'l'); % contains r-total, avoid unexpected energy contribution
model = changeRxnBounds(model,'EX18126',-0.001,'l'); % contains r-total, avoid unexpected energy contribution
model = changeRxnBounds(model,'SNK0069',-0.01,'l');
model = changeRxnBounds(model,'SNK0101',-0.01,'l');
model = changeRxnBounds(model,'SNK1002',-0.01,'l');
model = changeRxnBounds(model,'SNK1003',-0.01,'l');
model.rxnNotes(ismember(model.rxns,{'EX00089','EX00121','EX00132','EX00133','EX00185','EX00187','EX00208',...
                                    'EX00243','EX00257','EX00402','EX00486','EX00523','EX00535','EX00708',...
                                    'EX00777','EX01613','EX01801','EX02094','EX02415','EX05127','EX05402',...
                                    'EX05422','EX05697','EXC0152','EX01132','EX18125','EX18126','SNK0069',...
                                    'SNK0101','SNK1002','SNK1003'})) = ...
                                    {'Essential exchange reactions to ensure all reactions in the network can carry flux. For essential exchanges, we allow a default minimal flux of 0.01, except for EX18125 and EX18126 that imports R-total moieties. For EX18125 and EX18126, a lower flux limit of 0.001 was used.'};


% constraints to optimize the modeling of redox balance

% allow the following canonical nadph producing pathways unconstrained:PMID: 24805240
% RM00248: glu --> akg
% RC02736, RC01528: PPP
% RC00267, RM00267: idh-1, icit --> akg
% RM01220, RC01220: folate (MTHFD)
% RM00216: malate --> pyruvate

% To avoid loops, we assume directions for canonical pathways if there is
% a nadh version of them

% assume both the directions of MTHFD and MTHFD2 as reducing direction
% (making nadh or nadph) (PMID: 24805240)
model = changeRxnBounds(model,'RC01218',0,'l'); 
model = changeRxnBounds(model,'RC01220',0,'l');
model.rxnNotes(ismember(model.rxns,{'RC01218','RC01220'})) = {'Assuming both the directions of MTHFD and MTHFD2 as reducing direction (producing NADPH/NADH) in vivo, based on PMID: 24805240'};

% assuming reducing direction of GDH in vivo (PMID: 28208702)
model = changeRxnBounds(model,'RM00248',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model = changeRxnBounds(model,'RM00243',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
model.rxnNotes(ismember(model.rxns,{'RM00248','RM00243'})) = {'Assuming the direction of GDH as reducing direction (producing NADPH/NADH) in vivo, based on PMID: 28208702'};

% the following were constrianed to small flux (all other nadph
% producing reactions)
% of note, in most cases, the capacity of converting major reactant is
% not compromised by the constraint because there is another nadh
% version to support the flux
% RM00112 nnt-1 reaction was already constrained in above
model = changeRxnBounds(model,'RC04360',-0.01,'l');
model = changeRxnBounds(model,'RC04571',-0.01,'l');
model = changeRxnBounds(model,'RC07759',-0.01,'l');
model = changeRxnBounds(model,'RC01787',-0.01,'l');
model = changeRxnBounds(model,'RC05692',-0.01,'l');
model = changeRxnBounds(model,'RC01095',-0.01,'l');
model = changeRxnBounds(model,'RC02577',0.01,'u');
model = changeRxnBounds(model,'RC01041',-0.01,'l');
model = changeRxnBounds(model,'RM02566',0.01,'u');
model = changeRxnBounds(model,'RC00711',0.01,'u');
model = changeRxnBounds(model,'RM00711',0.01,'u');
model = changeRxnBounds(model,'RM00716',-0.01,'l');
model = changeRxnBounds(model,'RM03103',0.01,'u');
model = changeRxnBounds(model,'RC05623',-0.01,'l');
model = changeRxnBounds(model,'RC07140',0.01,'u');
model = changeRxnBounds(model,'RC08539',-0.01,'l');
model = changeRxnBounds(model,'RC00939',-0.01,'l');
model = changeRxnBounds(model,'RC02236',-0.01,'l');
model = changeRxnBounds(model,'RC01224',-0.01,'l');
model = changeRxnBounds(model,'RC00941',0.01,'u');
model = changeRxnBounds(model,'RC01904',-0.01,'l'); 
model = changeRxnBounds(model,'RC01431',-0.01,'l'); 
model = changeRxnBounds(model,'RC01759',-0.01,'l');
model = changeRxnBounds(model,'RC01481',-0.01,'l');
model = changeRxnBounds(model,'RM08759',0.01,'u'); 
model = changeRxnBounds(model,'RC08759',0.01,'u');
model = changeRxnBounds(model,'RM00706',0.01,'u');
model = changeRxnBounds(model,'RC00978',0.01,'u');
model = changeRxnBounds(model,'RC01415',0.01,'u');
model = changeRxnBounds(model,'RC08379',0.01,'u');
model = changeRxnBounds(model,'RC08383',0.01,'u');
model = changeRxnBounds(model,'RC03596',0.01,'u');
model = changeRxnBounds(model,'RC04940',-0.01,'l');
model = changeRxnBounds(model,'RC02082',-0.01,'l');
model = changeRxnBounds(model,'RC02697',0.01,'u');
model = changeRxnBounds(model,'RC03302',0.01,'u');
model = changeRxnBounds(model,'RC10059',0.01,'u');
model = changeRxnBounds(model,'RM03293',0.01,'u');
model.rxnNotes(ismember(model.rxns,{'RC04360','RC04571','RC07759','RC01787','RC05692','RC01095','RC02577',...
                                    'RC01041','RM02566','RC00711','RM00711','RM00716','RM03103','RC05623',...
                                    'RC07140','RC08539','RC00939','RC02236','RC01224','RC00941','RC01904',...
                                    'RC01431','RC01759','RC01481','RM08759','RC08759','RM00706','RC00978',...
                                    'RC01415','RC08379','RC08383','RC03596','RC04940','RC02082','RC02697',...
                                    'RC03302','RC10059','RM03293'})) = ...
                                    {'Restricting non-canonical NADPH production flux. Canonical reactions were selected based on PMID: 24805240 (incl. RM00248, RC02736, RC01528, RC00267, RM00267, RM01220, RC01220 and RM00216).'};


% other changes
model = changeRxnBounds(model,'RCC0139',0,'l'); % although BRENDA supports reversible, a sig. reverse flux is not likely feasible and this is the setting in human model
model.rxnNotes(strcmp(model.rxns,'RCC0139')) = {'Reaction RCC0139 was constrained to only allow forward flux (bounds = [0,1000]). This is because we found the reverse flux of this reaction can convert diphosphate (ppi) back to GTP without energy cost, forming a thermodynamically infeasible loop that recycles ppi produced in AA-tRNA synthesis to fuel GTP for protein synthesis without using energy (ATP). The modification of reversibility is supported by the reversibility annotation of the corresponding EC family (2.7.7.13) in the human model Human1.'};

model = changeRxnBounds(model,'TCM1071',0,'b'); 
model.rxnNotes(strcmp(model.rxns,'TCM1071')) = {'Block TCM1071 that is a misannotation in the model and should be not existing.'};


% add drains to uncouple biomass binding 
allPrecursors = model.mets(model.S(:, strcmp(model.rxns, 'BIO0107'))<0);
for i = 1:length(allPrecursors)
    model = addReaction(model,['DMN_',allPrecursors{i}],'reactionFormula',[allPrecursors{i},' ->'],'geneRule', 'NA','printLevel',1);
    model.rxnNotes(end) = {'New reaction added to uncouple the draining fluxes of different biomass precursors.'};
end


% add redox demand obj
model = addReaction(model,'OBJ_nadph','reactionFormula','nadph[c] -> nadp[c]','geneRule', 'NA','printLevel',1);
model.rxnNotes(end) = {'New reaction added to use as objective function of the redox production (not included in the CR model - keep for reproducibility purpose only).'};

% total lipid objectives
model = addReaction(model,'DMN_pail45p','reactionFormula','pail45p[c] ->','geneRule', 'NA','printLevel',1);
model.rxnNotes(end) = {'New reaction added to use as objective function of the pail45p production.'};

% total energy objectives
model = addReaction(model,'DMN_ISclustprot','reactionFormula','ISclustprot[m] -> apoISclustprot[m]','geneRule', 'NA','printLevel',1);
model.rxnNotes(end) = {'New reaction added to use as objective function of the iron-sulfur cluster production.'};
model = addReaction(model,'DMN_focytC','reactionFormula','focytC[m] -> apocytc[m]','geneRule', 'NA','printLevel',1);
model.rxnNotes(end) = {'New reaction added to use as objective function of the cytochromeC production.'};

%% save
model_out = struct();
model_out.rxnID  = model.rxns;
model_out.rxnFormula  = printRxnFormula(model,model.rxns,false);
model_out.LB  = model.lb;
model_out.UB  = model.ub;
model_out.note  = model.rxnNotes;
model_out = struct2table(model_out);

writetable(model_out,'output/iCEL1314_constriants.csv');
