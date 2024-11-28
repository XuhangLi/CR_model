function [model] = addDefaultConstraint(model,VitLim, essenlipidLim)
    % glucose limiting cannot be achived - the formal analysis may need to
    % use the real exchange flux contraint from NCI-60 to reasonably setup
    % the metabolism

    % set up the constraints for perturb-seq culture media 
    % K562 cells were grown in RPMI-1640 with 25 mM HEPES, 2.0 g/l NaHCO3, 
    % and 0.3 g/l L-glutamine supplemented with 10% FBS, 2 mM glutamine, 
    % 100 units/ml penicillin, and 100 μg/ml streptomycin. 
    % https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html

    % close all exchange 
    model.lb(findExcRxns(model)) = 0;

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
            else % add sink
                model = addSinkReactions(model,[myrxn,'[c]'],flux_lb(i, strcmp(cellLabel, 'LINE_K562')),1000);
            end
        end
    end
    
    % specially, the following trace nutrients are limiting biomass
    % production and may reflect inaccurate biomass composition; we allow
    % relaxed uptake of them to avoid such limits 
    model = changeRxnBounds(model,'MAR09145',VitLim,'l'); % pantothenate   

  
    % fill in default bounds for metabolites that is not measured but
    % important (e.g., vitamin)

    % the exchange reactions that we have real data for is commented out

    % first allow trace amount of essential metabolites (but allow it for
    % biomass/accounting for FBS, vitamin, etc)
    % these are known to be in the media
    
    % model = changeRxnBounds(model,'MAR09083',VitLim,'l'); % choline 
    % model = changeRxnBounds(model,'MAR09378',VitLim,'l'); % nictinomide  
    
    model = changeRxnBounds(model,'MAR09361',VitLim,'l'); % inositol 
    model = changeRxnBounds(model,'MAR09146',VitLim,'l'); % folate 
    model = changeRxnBounds(model,'MAR09144',VitLim,'l'); % pyridoxine  
    model = changeRxnBounds(model,'MAR09159',VitLim,'l'); % thiamine  
    model = changeRxnBounds(model,'MAR09109',VitLim,'l'); % biotin   
    model = changeRxnBounds(model,'MAR09143',VitLim,'l'); % riboflavin   
    model = changeRxnBounds(model,'MAR09269',VitLim,'l'); % B12   
    % all other cofactors in cofactor pool avaiable in very trace amout 
    cofactors = model.metNames(model.S(:,strcmp(model.rxns,'MAR10065')) < 0);
    cofactorEx = regexprep(cofactors, ' \[.+\]$',' [Extracellular]');
    cofactorExRxns = intersect(model.rxns(any(model.S(ismember(model.metNames,cofactorEx),:) < 0)), model.rxns(findExcRxns(model)));
    % skip ones with existing bonds (from data)
    cofactorExRxns = setdiff(cofactorExRxns, model.rxns(model.lb ~= 0));
    model = changeRxnBounds(model,cofactorExRxns,VitLim,'l');   
    % vitamins that is transformed (so not in cofactorEx) but needed 
    model = changeRxnBounds(model,{'MAR09151','MAR09154'},VitLim,'l'); % vitE
    model = changeRxnBounds(model,'MAR10492',VitLim,'l'); % retinol   


    % second, allow essential FA
    % they are C18 (three times carbon for glucose. To avoid significant carbon
    % flux (<1% glucose) ~1.5/3/100 = 0.005;
    model = changeRxnBounds(model,'MAR09035',essenlipidLim,'l'); % linoneate
    model = changeRxnBounds(model,'MAR09036',essenlipidLim,'l'); % linolenate


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


    % next a few aa are not measured, use a number from similar AA
    model = changeRxnBounds(model,'MAR09038',-0.01,'l'); % histidine; % a conservative number based on the flux of other essential AA (~0.01-0.03)
    model = changeRxnBounds(model,'MAR09363',-0.01,'l'); % cystine 

    % glutathione is also present in the medium. Let's use the AA minimal
    % flux as a proxy 
    model = changeRxnBounds(model,'MAR09351',-0.01,'l'); % glutathione 

    % set obj
    model = changeObjective(model,'MAR13082');% we stay with the default biomass instead of total biomass (MAR10024) for now



end
