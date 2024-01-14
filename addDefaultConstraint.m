function [model] = addDefaultConstraint(model,type)
%%
if strcmp(type,'default')
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',0,'u');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
    model = changeObjective(model,'BIO0010');
elseif strcmp(type,'nutritionalFree@1000')
    model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = -1000;
    model = changeObjective(model,'BIO0010');
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',0,'u');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1000,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
elseif strcmp(type,'nutritionalFree@1')
    model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = -1;
    model = changeObjective(model,'BIO0010');
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',0,'u');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');
    model = changeRxnBounds(model,'RCC0005',1,'l');
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0,'u');
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');
    model = changeRxnBounds(model,'EXC0050',-1,'l');
    model = changeRxnBounds(model,'EXC0050',0,'u');
elseif strcmp(type,'minimalExchange@1')
    model.lb(cellfun(@(x) strcmp(x(1:2),'EX'),model.rxns)) = 0;
    model.lb(cellfun(@(x) strcmp(x(1:3),'SNK'),model.rxns)) = 0;
    model = changeObjective(model,'BIO0010');
    model = changeRxnBounds(model,'DMN0033',0,'l');
    model = changeRxnBounds(model,'DMN0033',0.01,'u');
    model = changeRxnBounds(model,'RCC0005',10,'l'); % assume NGAM as 10x bacterial uptake
    model = changeRxnBounds(model,'RCC0005',1000,'u');
    model = changeRxnBounds(model,'RM00112',0,'l');
    model = changeRxnBounds(model,'RM00112',0.01,'u'); % nnt-1 reaction; assume producing nadph and assume not massive production (since H+ gradient coupling was not reconstructed)
    model = changeRxnBounds(model,'RM04432',-1000,'l');
    model = changeRxnBounds(model,'RM04432',1000,'u');
    model = changeRxnBounds(model,'RMC0005',0,'u');
    model = changeRxnBounds(model,'RMC0005',0,'l');

    % Essential exchanges:
    model = changeRxnBounds(model,'EXC0050',-1,'l');
    model = changeRxnBounds(model,'EX00001',-1000,'l');
    model = changeRxnBounds(model,'EX00007',-1000,'l');
    model = changeRxnBounds(model,'EX00009',-1000,'l');
    model = changeRxnBounds(model,'EX00080',-1000,'l');

    % allowed storage (to simulate genes related to storage ultilization)
    model = changeRxnBounds(model,'SNK0012',-0.01,'l');
    model = changeRxnBounds(model,'SNK0013',-0.01,'l');
    model = changeRxnBounds(model,'SNK0014',-0.01,'l'); % start with 0.01 TAG and this can be lower if it causes irrealisitic carbon influx

 
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
    
    % assuming reducing direction of GDH in vivo (PMID: 28208702)
    model = changeRxnBounds(model,'RM00248',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study
    model = changeRxnBounds(model,'RM00243',0,'l'); % assume glutamate dehydrogenase only in forward direction (making nh4) in our study


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

    % block the thermodynamically infeasible loops that can use the low-energy
    % bound in ATP as a high-energy bound
    % the massive PPI produced in tRNA synthesis can be converted to GTP in a
    % loop that further support protein synthesis; this is an infeasible loop
    % to bypass real energy demand; we prevent this loop
    model = changeRxnBounds(model,'RCC0139',0,'l'); % although BRENDA supports reversible, a sig. reverse flux is not likely feasible and this is the setting in human model

    
end
end
