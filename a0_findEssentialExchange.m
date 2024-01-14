initCobraToolbox

%% 
% we will generate a minimal exchange model that only essential exchange
% reactions are allowed (including sinks for storage molecules)


%%
load('iCEL1314.mat');
worm = addDefaultConstraint(model,'nutritionalFree@1');
% worm = changeRxnBounds(worm,'RCC0005',10,'l');
% additional SNK for nutritional free state
worm = changeRxnBounds(worm,'SNK0069',-1,'l');
worm = changeRxnBounds(worm,'SNK0101',-1,'l');
worm = changeRxnBounds(worm,'SNK1002',-1,'l');
worm = changeRxnBounds(worm,'SNK1003',-1,'l');
% make the epsilon vector
[epsilon_f,epsilon_r, capacity_f, capacity_r] = makeEpsilonSeq(worm,worm.rxns,0.01,0.5);
%parsedGPR = GPRparser_xl(worm);% Extracting GPR data from model
%worm.parsedGPR = parsedGPR;

%% 1. check if all non-exchange rxns can be opened 
worm = addDefaultConstraint(model,'nutritionalFree@1000');
% worm = changeRxnBounds(worm,'RCC0005',10,'l');
% additional SNK for nutritional free state
worm = changeRxnBounds(worm,'SNK0069',-1000,'l');
worm = changeRxnBounds(worm,'SNK0101',-1000,'l');
worm = changeRxnBounds(worm,'SNK1002',-1000,'l');
worm = changeRxnBounds(worm,'SNK1003',-1000,'l');

inds = cellfun(@(x) contains(x,'EX'),model.rxns) | cellfun(@(x) contains(x,'DM'),model.rxns) | cellfun(@(x) contains(x,'SNK'),model.rxns) | cellfun(@(x) contains(x,'BIO'),model.rxns) | cellfun(@(x) contains(x,'UP'),model.rxns);
%add exceptions for nonGPR transporters; this is to avoid that some
%exchange becoming essential only because of the nonGPR transporter
exceptions = {'TCE0766','TCE5359','TCM0170','TCE0739','TCE1158','TCE0558','TCE0374','TCE0591','TCE5074','TCE5073','TCE5109'}
inds(ismember(model.rxns,exceptions)) = true;
mimic = 3*ones(length(model.rxns),1);
mimic(inds) = 2;

[solution, MILProblem,Hreactions,Lreactions] = iMAT_xl(worm, mimic, epsilon_f, epsilon_r,1.1,2.9);
flux = solution.full(1:length(model.rxns));
onRxns = flux >= epsilon_f-1e-8 | flux <= -epsilon_r+1e-8;
failOpen(:,1) = setdiff(Hreactions,model.rxns(onRxns));
[A B] = ismember(failOpen(:,1),model.rxns);
capaci_f = (capacity_f(B(A)));
capaci_r = (capacity_r(B(A)));
failOpen(:,2) = mat2cell(capaci_f,ones(length(failOpen(:,1)),1),1);
failOpen(:,3) = mat2cell(capaci_r,ones(length(failOpen(:,1)),1),1);
%% find the minimal Exchange
MILProblem2 = solution2constraint(MILProblem,solution);
RHindex = find(cellfun(@(x) contains(x,'EX'),model.rxns)); %>= 0
Hreactions = model.rxns(RHindex);
epsilon_f_sorted = zeros(length(RHindex));
S = MILProblem2.A;
lb = MILProblem2.lb;
ub = MILProblem2.ub;
% Creating A matrix
A = sparse(size(S,1)+length(RHindex),size(S,2)+length(RHindex));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end

for i = 1:length(RHindex)
    A(i+size(S,1),RHindex(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon_f_sorted(i);
end

% Creating csense
csense1(1:size(S,1),1) = MILProblem2.csense;
csense2(1:length(RHindex),1) = 'G';
csense = [csense1;csense2];

% Creating lb and ub
lb_y = zeros(length(RHindex),1);
ub_y = ones(length(RHindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];

% Creating c
c_v = zeros(size(S,2),1);
c_y = ones(length(RHindex),1);
c = [c_v;c_y];

% Creating b
b_s = MILProblem2.b;
lb_rh = lb(RHindex);
b = [b_s;lb_rh];

% Creating vartype
vartype1(1:size(S,2),1) = MILProblem2.vartype;
vartype2(1:length(RHindex),1) = 'B';
vartype = [vartype1;vartype2];

MILPproblem2.A = A;
MILPproblem2.b = b;
MILPproblem2.c = c;
MILPproblem2.lb = lb;
MILPproblem2.ub = ub;
MILPproblem2.csense = csense;
MILPproblem2.vartype = vartype;
MILPproblem2.osense = -1;
MILPproblem2.x0 = [];

solution2 = solveCobraMILP_XL(MILPproblem2, 'timeLimit', 7200, 'printLevel', 1);
%% get minimal EX
EXrxns = model.rxns(RHindex);
closedEX = EXrxns(solution2.full(RHindex) >= -1e-8);
minimalEX = setdiff(EXrxns,closedEX)

% Essential exchanges:
% worm = changeRxnBounds(worm,'EXC0050',-1000,'l');
% worm = changeRxnBounds(worm,'EX00001',-1000,'l');
% worm = changeRxnBounds(worm,'EX00007',-1000,'l');
% worm = changeRxnBounds(worm,'EX00009',-1000,'l');
% worm = changeRxnBounds(worm,'EX00080',-1000,'l');
% 
% Additional exchanges:
% worm = changeRxnBounds(worm,'EX00089',-1000,'l');
% worm = changeRxnBounds(worm,'EX00121',-1000,'l');
% worm = changeRxnBounds(worm,'EX00132',-1000,'l');
% worm = changeRxnBounds(worm,'EX00133',-1000,'l');
% worm = changeRxnBounds(worm,'EX00185',-1000,'l');
% worm = changeRxnBounds(worm,'EX00187',-1000,'l');
% worm = changeRxnBounds(worm,'EX00208',-1000,'l');
% worm = changeRxnBounds(worm,'EX00243',-1000,'l');
% worm = changeRxnBounds(worm,'EX00257',-1000,'l');
% worm = changeRxnBounds(worm,'EX00402',-1000,'l');
% worm = changeRxnBounds(worm,'EX00486',-1000,'l');
% worm = changeRxnBounds(worm,'EX00523',-1000,'l');
% worm = changeRxnBounds(worm,'EX00535',-1000,'l');
% worm = changeRxnBounds(worm,'EX00708',-1000,'l');
% worm = changeRxnBounds(worm,'EX00777',-1000,'l');
% worm = changeRxnBounds(worm,'EX01613',-1000,'l');
% worm = changeRxnBounds(worm,'EX01801',-1000,'l');
% worm = changeRxnBounds(worm,'EX02094',-1000,'l');
% worm = changeRxnBounds(worm,'EX02415',-1000,'l');
% worm = changeRxnBounds(worm,'EX05127',-1000,'l');
% worm = changeRxnBounds(worm,'EX05402',-1000,'l');
% worm = changeRxnBounds(worm,'EX05422',-1000,'l');
% worm = changeRxnBounds(worm,'EX05697',-1000,'l');
% worm = changeRxnBounds(worm,'EXC0152',-1000,'l');
% worm = changeRxnBounds(worm,'EX01132',-1000,'l');
% worm = changeRxnBounds(worm,'EX18125',-1000,'l');
% worm = changeRxnBounds(worm,'EX18126',-1000,'l');
% worm = changeRxnBounds(worm,'SNK0069',-1000,'l');
% worm = changeRxnBounds(worm,'SNK0101',-1000,'l');
% worm = changeRxnBounds(worm,'SNK1002',-1000,'l');
% worm = changeRxnBounds(worm,'SNK1003',-1000,'l');
