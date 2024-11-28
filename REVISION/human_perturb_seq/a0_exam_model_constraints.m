load('input/human_model/Human-GEM_cobra.mat'); % we found ihuman doesnt have csense field - so we use cobra 

model = addDefaultConstraint(model,-0.001, -0.005); % arbiturary flux <1% glucose carbon
% has to figure out individual essential FA
sol = optimizeCbModel(model)

% growth rate close to exp. data (0.035)

%% check limiting nutrient 

allcons = model.rxns(model.lb ~= 0 & model.lb ~= -1000);
lims = [];
for i = 1:length(allcons)
    test = changeRxnBounds(model,allcons{i}, -10,'l');
    sol = optimizeCbModel(test);
    if sol.f > 0.047
        lims = [lims; allcons(i)];
        lims
    end

end
% tryptophan is now limiting (exp flux - good).
%%
model.lb(model.c==1) = sol.f;
flux = minimizeModelFlux_XL(model);
listRxn(model,flux,'acetyl-CoA [Mitochondria]')
% PDH is also active 

%%
listRxn(model,flux,'pyruvate [Cytosol]')


%%
model = changeObjective(model,'MAR03964');
sol = optimizeCbModel(model)

%%
seq = 0:5:50;
biomass = [];
for i = 1:length(seq)
    model = changeRxnBounds(model,'MAR09034',-seq(i),'l'); % FA pool for some uncommon essential FA
    sol = optimizeCbModel(model);
    biomass(i) = sol.f;
end
figure
plot(seq, biomass)

%%
model.metNames(sol.y > 0)
sol
%%

a = model.metNames(sol.y > 10);
a(contains(a,'Extracellular'))
intersect(a, allFA)
%%
model = addReaction(model,'test','reactionFormula',['MAM01589c[c]',' ->'],'geneRule', 'NA','printLevel',1);
model = changeObjective(model,'test');
sol = optimizeCbModel(model)
a = model.metNames(sol.y > 10);
a(contains(a,'Extracellular'));
%%
seq = [0:0.001:0.01,0.02:0.01:0.1,0.2:0.1:1];
biomass = [];
for i = 1:length(seq)
    model = changeRxnBounds(model,'MAR13039',-seq(i),'l'); % FA pool for some uncommon essential FA
    sol = optimizeCbModel(model);
    biomass(i) = sol.f;
end
%% glucose
seq = 10:10:30;
biomass = [];
for i = 1:length(seq)
    model = changeRxnBounds(model,'MAR09034',-seq(i),'l'); % FA pool for some uncommon essential FA
    sol = optimizeCbModel(model);
    biomass(i) = sol.f;
end
figure
plot(seq, biomass)
%% vit
seq = [0:0.002:0.05, 0.02:0.02:0.1, 0.2:0.2:1, 2:2:10];
biomass = [];
for i = 1:length(seq)
    model = addDefaultConstraint(model,-seq(i), -1); % vitamin not limiting 
    sol = optimizeCbModel(model)   
    biomass(i) = sol.f;
end
figure
plot(seq, biomass)

%% FA
seq = [0:0.002:0.05, 0.02:0.02:0.1, 0.2:0.2:1, 2:2:10];
biomass = [];
for i = 1:length(seq)
    model = addDefaultConstraint(model,-0.1, -seq(i)); % vitamin not limiting 
    sol = optimizeCbModel(model)   
    biomass(i) = sol.f;
end
figure
plot(seq, biomass)

