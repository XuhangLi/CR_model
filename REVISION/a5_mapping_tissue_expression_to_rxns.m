% we used the GPR parsing methods from FPA/eFPA

% calculate the penalty matrix given expression data. It translates the
% expression level of genes to level of reactions by a set of GPR parsing
% rules. In brief, we taking the "sum" of all "OR" gated genes in a single
% GPR and "min" of "AND" ones. Complex GPRs are handled specially by
% converting the GPR annotation to functionally "AND" gated blocks, where
% all "OR" gated genes are added up. In the convertion to blocks, if a
% "AND" gate is nested in an "OR", the minimal value will be taken. For
% detailed explanation, see supplemental information. 
addpath scripts\
%% load gene data
tissueTPM = readtable('./../../2_DE/clustering_python/scRNA_data_processing/tissueTPM_combined.csv');
IDtbl = readtable('input\IDtbl.csv');
tissueTPM = tissueTPM(ismember(tissueTPM.Var1, IDtbl.WormBase_Gene_ID),:);
[A B] = ismember(tissueTPM.Var1, IDtbl.WormBase_Gene_ID);
tissueTPM.Var1 = IDtbl.ICELgene(B(A));
load('./../CR_model_final_run/input/iCEL1314.mat');

%% step1: mapping the gene-centric expression levels to reactions
levels = cell(length(model.rxns),size(tissueTPM, 2)-1);
status = -1*ones(length(model.rxns),size(tissueTPM, 2)-1);
fprintf('mapping expression levels to penalties...\n');
fprintf(['\n' repmat('.',1,size(tissueTPM, 2)-1) '\n\n']);
for i = 2:size(tissueTPM, 2)
    expression = tissueTPM{:,i} + 1; % add one pseudo count for stability
    % expression.value = expression.value + 1; % we dont add pseudocount - it will be preprocessing if needed!
    [levels(:,i-1), status(:,i-1)] = gene_to_reaction_levels(model, tissueTPM.Var1, expression, @min, @(x,y)(x+y));%GPR parser to convert the expression of genes to levels of "AND" gated blocks
    fprintf('\b|\n');%for simple progress monitor
end

%% option 1: absolute expression 
% since we want to first look at the absolute expression (to avoid biased
% interpretation), we will take the minimal of all functional "AND"
% directly, instead of do a "block-wise normalization" as in FPA 
merged_levels = [];
for i = 1:size(levels)
    for j = 1:size(levels, 2)
       merged_levels(i,j) = min(levels{i,j});
    end
end

% save 
merged_levels = array2table(merged_levels);
merged_levels.Properties.VariableNames = tissueTPM.Properties.VariableNames(2:end);
merged_levels.Properties.RowNames = model.rxns;
% this tpm has pseudo count added
writetable(merged_levels, 'input/tissue_TPM_parsed.csv','WriteRowNames',true);

%% option 2: relative expression 
normalizedLevel = nan(size(levels,1),size(levels,2));%the levels normalized to super condition
penalty = ones(size(levels,1),size(levels,2));%set an default penalty as 1
for i = 1:length(model.rxns)
    % this is for error handling purposes. If some very special GPR format
    % appears and the parser fails to correctly handle it, the number of
    % blocks for the same reaction in different conditions is likely
    % different. So, we check the consistance first to make sure the GPR
    % parsing is as expected.
    lenL = zeros(size(levels, 2),1);
    for j=1:size(levels, 2)
        lenL(j)=length(levels{i,j});
    end
    if any(lenL ~= lenL(1)) %some GPR parsing error
        error('GPR length unequal');
    else
        % normalize every block
        stackM = nan(lenL(1),size(levels, 2)); %stack the levels of blocks in a matrix (it was in cell variable originally)
        for z = 1:lenL(1)
            for s = 1:size(levels, 2)
                stackM(z,s) = levels{i,s}(z);
            end
        end
        % first normalize the matrix to the super condition (highest value)
        stackM = stackM ./ max(stackM,[],2);
        % then, pick the minimal value of all blocks (because they are "AND" connected
        normalizedLevel(i,:) = min(stackM,[],1);
        if ~any(isnan(normalizedLevel(i,:))) %if no NaN occurs
            penalty(i,:) = 1./normalizedLevel(i,:);
        % if all the levels (in different conditions) are NaN, this
        % reaction must have no data. We will keep the default penalty (which is 1) for
        % it 
        elseif any(isnan(normalizedLevel(i,:))) && ~all(isnan(normalizedLevel(i,:))) %if only some condition gives NaN
            error('partial None Penalty for some conditions, check!');
        end
    end
end

% save 
merged_levels = array2table(normalizedLevel);
merged_levels.Properties.VariableNames = tissueTPM.Properties.VariableNames(2:end);
merged_levels.Properties.RowNames = model.rxns;
writetable(merged_levels, 'input/tissue_relative_expression_of_reactions.csv','WriteRowNames',true);

