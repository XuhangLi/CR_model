function mytbl = listRxn4gene(model,flux,myGene)
% This is a simple flux tracker to show the flux around a metabolite
% according to a flux distribution
%
% USAGE:
%
%    mytbl = listRxn(model,flux,myMet)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    flux:              the flux distribution to inspect
%    myMet:             the metabolite of interest to list flux around
%
% OUTPUT:
%   mytbl:             	a table showing all the flux that produces or
%                       comsumes myMet
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020
myrxns = model.rxns(model.rxnGeneMat(:,strcmp(model.genes,myGene))==1);
totalflux = 0;
for i = 1:length(myrxns)
    mytbl(i,1) = myrxns(i);
    mytbl(i,2) = {flux(strcmp(model.rxns,myrxns{i}))};
    mytbl(i,3) = printRxnFormula(model, myrxns{i},0);
    totalflux = totalflux + abs(flux(strcmp(model.rxns,myrxns{i})));
end
mytbl = sortrows(mytbl,3);
fprintf('total flux is %3f\n',totalflux);
end