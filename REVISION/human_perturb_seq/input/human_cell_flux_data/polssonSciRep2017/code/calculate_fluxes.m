function [fluxes, modelCore,lineLabels, modelsAll] = calculate_fluxes()

%% NOTES:

%This function accompanies the manuscript "Systems biology analysis of 
%drivers underlying hallmarks of cancer cell metabolism"
%Zielinski et. al. 

% %Required third party files and code are included in the folder.
%Required files and code:
%File: varSet1.mat
%File: modelCore.mat
%Code: optimizeCbModel_Ratio.m
%Code: solveCobraQP_custom.m
%Code: COBRA toolbox  https://opencobra.github.io/
%Solver: Gurobi 5 http://www.gurobi.com/

%LP and QP solver is required. The code was tested with Gurobi 5, which is free
%for academic use (http://www.gurobi.com/). Other solvers may work but the
%changeCobraSolver lines below must be altered accordingly.


%% Start

%Initialize solvers: see changeCobraSolver in the COBRA toolbox for useable
%solvers
changeCobraSolver('gurobi5','LP')
changeCobraSolver('gurobi5','QP')

scriptName = mfilename('fullpath');
[currentpath, filename, fileextension]= fileparts(scriptName);

%NOTE: IF THIS LINE THROWS AN ERROR, YOU CAN MANUALLY INSERT THE
%APPROPRIATE ABSOLUTE PATH FOR THE CORRESPONDING VARIABLES
varName1 = [currentpath,'./','varSet1.mat'];
varName2 = [currentpath,'./','modelCore.mat'];


%% Variable import

%Import data (exchanges, growth rates)
load(varName1);
% chr;
dataDoubled; %In fmol/cell/hr
grMaxes; %in /hr
grMins; %in /hr
lineLabels;
linesDoubled;
masses; % g/cell
modelRecon2;
o2Max; %in mmol/gDW/hr
o2Min; %in mmol/gDW/hr
rxnsData;
proteinCorrected; %g/cell
proteinMass; %g/cell
proteinSecreted; %g/cell
volumes; %in pL

%Import the core cancer model
load(varName2 )
% modelCore = addReaction(modelCore,'DM_nadh[c]',{'nadh[c]','nad[c]','h[c]'},[-1 1 1],0,0,1000,1);

%% Process data
%Remove data for SR and MDA-N from all variables, since data is either
%missing or unreliable from these lines
linesToRemove = {'LINE_MDA_MB_468'};
indsLabels = cell2mat(cellfun(@(x) find(strcmp(x,lineLabels)),linesToRemove,'UniformOutput',false));
linesToRemoveAlt = {'MDA-MB-468'};
indsDoubled = cell2mat(cellfun(@(x) find(strcmp(x,linesDoubled)),linesToRemoveAlt,'UniformOutput',false));

chr(indsLabels) = [];
dataDoubled(:,indsDoubled) = [];
grMaxes(indsLabels) = [];
grMins(indsLabels) = [];
lineLabels(indsLabels) = [];
linesDoubled(indsDoubled) = [];
masses(indsLabels) = [];
o2Max(indsLabels) = [];
o2Min(indsLabels) = [];
proteinMass(indsLabels) = [];
proteinSecreted(indsLabels) = [];
volumes(indsLabels) = [];

dataDoubledAdj = dataDoubled;

%Calculate the bounds
dataMean = [];
dataStd = [];
for i = 1:2:size(dataDoubledAdj,2)
    dataMeanCur = mean([dataDoubledAdj(:,i),dataDoubledAdj(:,i+1)],2);
    dataStdCur = std([dataDoubledAdj(:,i),dataDoubledAdj(:,i+1)],0,2);
    dataMean = [dataMean,dataMeanCur];
    dataStd = [dataStd,dataStdCur];
end

dataLow = dataMean-2*dataStd;
dataHigh = dataMean+2*dataStd;


%% Calculate flux states for minimal models

massFracProt = 0.70;
massesCorrected = proteinCorrected/massFracProt; 

%Make the unconstrained bounds more reasonable since scale is smaller in
%mass units
modelCore.lb(modelCore.lb<-1000) = -1000;
modelCore.ub(modelCore.ub>1000) = 1000;

rxnsOpen = {'EX_o2(e)';'EX_cl(e)';'EX_h2o(e)';'EX_k(e)';'EX_na1(e)';'EX_pi(e)';...
    'EX_chol(e)'};

rxnsConstrained = {'EX_glc(e)';'EX_lac_L(e)';'EX_gln_L(e)';'EX_glu_L(e)';'EX_asp_L(e)';...
    'EX_asn_L(e)';'EX_pro_L(e)';'EX_arg_L(e)';'EX_ala_L(e)';'EX_ser_L(e)';'EX_gly(e)';...
    'EX_lys_L(e)';'EX_trp_L(e)';'EX_leu_L(e)';'EX_tyr_L(e)';'EX_phe_L(e)';...
    'EX_ile_L(e)';'EX_val_L(e)';'EX_thr_L(e)';'DM_gudac_c_';...
    'EX_orn(e)';'EX_cit(e)';'EX_mal_L(e)'};

%Reset the constrained reactions
senseConstr = '';
matConstr = [];

%TESTING ROS PRODUCTION AS O2 
% %Set up equivalent O2 demand constraints
% ratioROS = 0.01;
% %Mitochondrial ROS production
% vecConstrMit = zeros(1,length(modelCore.rxns));
% vecConstrMit(find(strcmp('DM_o2[m]',modelCore.rxns))) = 1;
% rxnsNadphMit = {'NADH2_u10m';'CYOR_u10m';'CYOOm2';'SUCD1m';'ICDHxm';'ICDHyrm';'PDHm';'AKGDm'};
% signNadphMit = [1;1;1;1;1;1;1;1];
% for i = 1:length(rxnsNadphMit)
%     vecConstrMit(find(strcmp(rxnsNadphMit{i},modelCore.rxns))) = -signNadphMit(i)*ratioROS;
% end
% matConstr = [matConstr;vecConstrMit];
% senseConstr = [senseConstr;'L'];

%Add ICDHxm ICDHyrm even flux split constraint
vecConstrIcdh = zeros(1,length(modelCore.rxns));
vecConstrIcdh(find(strcmp('ICDHyrm',modelCore.rxns))) = 1;
vecConstrIcdh(find(strcmp('ICDHxm',modelCore.rxns))) = -1;
matConstr = [matConstr;vecConstrIcdh];
senseConstr = [senseConstr;'E'];

%Add GLUDxm GLUDym even flux split constraint
vecConstrGlud = zeros(1,length(modelCore.rxns));
vecConstrGlud(find(strcmp('GLUDym',modelCore.rxns))) = 1;
vecConstrGlud(find(strcmp('GLUDxm',modelCore.rxns))) = -1;
matConstr = [matConstr;vecConstrGlud];
senseConstr = [senseConstr;'E'];

%Enforce PPP 10% split, because otherwise things just don't work
vecConstrPPP = zeros(1,length(modelCore.rxns));
vecConstrPPP(find(strcmp('HEX1',modelCore.rxns))) = 0.10;
vecConstrPPP(find(strcmp('G6PDH2r',modelCore.rxns))) = -1;
matConstr = [matConstr;vecConstrPPP];
senseConstr = [senseConstr;'E'];

%Enforce PDHm/PCm 10% split, because otherwise things just don't work
vecConstrPDH = zeros(1,length(modelCore.rxns));
vecConstrPDH(find(strcmp('PDHm',modelCore.rxns))) = 0.10;
vecConstrPDH(find(strcmp('PCm',modelCore.rxns))) = -1;
matConstr = [matConstr;vecConstrPDH];
senseConstr = [senseConstr;'E'];


%Convert from fmol/cell/hr to mmol/gDW/hr using the cell masses
%First convert from fmol/cell/hr to mmol/cell/hr
dataLowMmol = dataLow/(10^12);
dataHighMmol = dataHigh/(10^12);
%Then convert to mmol/gDW/hr using the cell mass data
dataLowMass = bsxfun(@times,dataLowMmol,massesCorrected'.^-1);
dataHighMass = bsxfun(@times,dataHighMmol',massesCorrected.^-1)';
%Extract data for just the exchanges we want to use
rxnsConstrained;
dataLowMassRed = zeros(length(rxnsConstrained),size(dataLow,2));
dataHighMassRed = zeros(length(rxnsConstrained),size(dataHigh,2));
for i = 1:length(rxnsConstrained)
    curRxnInd = find(strcmp(rxnsConstrained{i},rxnsData));
    for j = 1:size(dataLow,2)
        dataLowMassRed(i,j) = dataLowMass(curRxnInd,j);
        dataHighMassRed(i,j) = dataHighMass(curRxnInd,j);
    end
end

%Set as a hard value
massFracProt = 0.70;

%Get mass fractions of each non DNA or protein component of biomass
massFracDNA = 0.049833754; %Back-calculated from diploid
componentsRed = {'rna';'glycogen';'resolved lipids';'resolved small molecules'};
componentFracsRed = [0.202445098;0.208094729;0.300681932;0.051111362];
components = ['protein';'dna';componentsRed];
componentFracs = zeros(length(lineLabels),length(components));
for i = 1:length(lineLabels)
    componentFracs(i,1) = massFracProt;
    componentFracs(i,2) = (1-massFracProt)*chr(i)/46*massFracDNA;
    for j = 1:length(componentsRed)
        componentFracs(i,2+j) = (1-massFracProt)*(1-componentFracs(i,2))*componentFracsRed(j);
    end
end

%Now calculate the mass-scaled biomass coefficients
%Calculate biomass value for each cell line in mmol/gDW cell
% Set up biomass reaction
biomassMets = {'trp_L[c]';'ala_L[c]';'arg_L[c]';'asn_L[c]';'asp_L[c]';'gln_L[c]';...
    'glu_L[c]';'gly[c]';'ile_L[c]';'leu_L[c]';'lys_L[c]';...
    'phe_L[c]';'pro_L[m]';'ser_L[c]';'thr_L[c]';'tyr_L[c]';'val_L[c]';...
    'datp[c]';'dctp[c]';'dgtp[c]';'dttp[c]';'ctp[c]';'gtp[c]';'utp[c]';'atp[c]';...
    'glygn2[c]';'sphmyln_hs[c]';'chsterol[r]';'xolest_hs[r]';'mag_hs[c]';'dag_hs[c]';...
    'pail_hs[c]';'pe_hs[c]';'ps_hs[c]';'pchol_hs[c]';'lpchol_hs[c]';'clpn_hs[c]';...
    'pa_hs[c]';'hdcea[c]';'hdca[c]';'ocdcea[c]';'ocdca[c]'};

biomassCompMassFracs = [0.0106;0.0597;0.0854;0.0527;0.0426;0.0584;0.0709;0.0411;...
    0.0517;0.0887;0.102;0.0359;0.033;0.048;0.0546;0.035;0.0546;0.304;...
    0.187;0.214;0.295;0.17;0.322;0.323;0.184;1;0.0533;0.127;0.116;0.00562;0.0069;...
    0.0477;0.149;0.059;0.198;0.013;0.0103;0.0846;0.00944;0.0238;0.0806;0.0158];

biomassCompInds = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;3;3;3;3;4;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5];

%molecular weights in g/mol
biomassMWs = [186.21; 71.08; 157.2; 114.1; 114.08; 128.13; 128.11; 57.05; 113.16;...
    113.16; 129.18; 147.18; 97.12; 87.08; 101.11; 163.18; 99.13; 312.2; 286.16;...
    328.2; 303.19; 304.18; 344.2; 305.16; 328.2; 1783.54; 463.62; 386.66; 413.67; 343;...
    594; 845; 727; 770; 769; 467.6; 1466.1; 706.2; 253; 255; 255; 283];

%ATP consumed in macromolecule synthesis for 1 gram of biomass, in mmol/gDW
maintenanceVals = [35;-35;-35;35;-35];

biomassVals = [];
for i = 1:length(lineLabels)
    curComps = componentFracs(i,:);
    curVals = zeros(1,length(biomassMets));
    for j = 1:length(biomassMets)
        curCompInd = biomassCompInds(j);
        curCompVal = curComps(curCompInd);
        curMW = biomassMWs(j);
        curVal = biomassCompMassFracs(j)*curCompVal/curMW*1000;
        curVals(j) = curVal;
    end
    curVals = [curVals,maintenanceVals'];
    
    biomassVals = [biomassVals;curVals];
end

maintenanceMets = {'atp[c]';'adp[c]';'pi[c]';'h2o[c]';'h[c]'};
biomassMets = [biomassMets; maintenanceMets];

%Set ATP demand at the real rate now
valATPM = 1.07;
modelCore = changeRxnBounds(modelCore,'DM_atp_c_',valATPM,'l');
modelCore = changeRxnBounds(modelCore,'DM_atp_c_',1000,'u');

%Make sure exchange inds are correct
rxnsConstrainedInds = cell2mat(cellfun(@(x) find(strcmp(x,modelCore.rxns)),rxnsConstrained,'UniformOutput',false));

%Now calculate fluxes
fluxes = [];
modelsAll = cell(length(lineLabels),1);
for i = 1:length(lineLabels)
%     i
%     curLine = lineLabels{i}
    curModel = modelCore;
    
    %Open oxygen, doesn't hit the bound
    curModel = changeRxnBounds(curModel,'EX_o2(e)',-20,'l');
    curModel = changeRxnBounds(curModel,'EX_o2(e)',0,'u');
    
    %Add in a biomass reaction
    curBiomassVals = biomassVals(i,:);
    curModel = addReaction(curModel,'biomass_NCI60',biomassMets,-curBiomassVals,0,0,1000,0);
    
    %Quick fix to make sure dimensions map - causes no problem, just a
    %little ugly
    if size(matConstr,2)~=size(curModel.S,2)
        matConstr = [matConstr,zeros(size(matConstr,1),1)];
    end
    
    %Set growth
    curGrowthMax = grMaxes(i);
    curGrowthMin = grMins(i);
    curModel = changeRxnBounds(curModel,'biomass_NCI60',curGrowthMax,'u');
    curModel = changeRxnBounds(curModel,'biomass_NCI60',curGrowthMin,'l');


    %Exchanges
    curExchangesLow = dataLowMassRed(:,i);
    curExchangesHigh = dataHighMassRed(:,i);
    for j = 1:length(rxnsConstrained)
        curExInd = rxnsConstrainedInds(j);
        if isempty(curExInd)
            rxnsConstrained(j)
        end
        curValHigh = curExchangesHigh(j);
        curModel.ub(curExInd) = curValHigh;
        curValLow = curExchangesLow(j);
        curModel.lb(curExInd) = curValLow;
    end
    
    %Fix reversibility information for the current model
    for j = 1:length(curModel.rxns)
        if curModel.lb(j)<0&&curModel.ub(j)>0
            curModel.rev(j) = 1;
        else
            curModel.rev(j) = 0;
        end
    end
    
    curModel = changeObjective(curModel,'DM_nadph[m]',1);
    curSol = optimizeCbModel_Ratio(curModel,'max',2,1,matConstr,senseConstr);

    fluxes = [fluxes,curSol.x];
    modelsAll{i} = curModel;
end
