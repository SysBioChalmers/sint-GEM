%repMets_analysis

%load model 
load('../../models/candida_intermedia/cintGEM.mat')
%Correct model grRules
model.grRules = strrep(model.grRules,'Candida_intermedia@','');
[grRules,rxnGeneMat] = standardizeGrRules(model,false);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;

cSources = {'gal' 'lac' 'cel' 'xyl'};
mkdir('../../results/reporter_metabolites')
multiDim_data = table();

%Exclude highly connected Mets
toExclude = {'carbon dioxide' 'ATP' 'ADP' 'AMP' 'phosphate' 'diphosphate' 'H+' ...
             'CMP' 'GDP' 'GTP' 'H2O' 'CTP' 'coenzyme A'};

for cSource=cSources
    %load data
    str = cSource{1};
    disp(str)
    disp(' ')
    DEresults = readtable(['../../results/RNA_DE_analysis/RNA_DE_glu_vs_' str '.txt'],'delimiter','\t');
    DEmapping = readtable(['../../results/RNA_DE_analysis/RNA_2_model_glu_vs_' str '.txt'],'delimiter','\t');
    nonZeros  = sum(DEmapping.counts_ref+DEmapping.counts_Csrc,2)>0;
    DEmapping = DEmapping(nonZeros,:);
    %
    if isempty(multiDim_data)
        multiDim_data.genes = DEmapping.modelGenes;
    end
    eval(['multiDim_data.pVal_' str '=DEmapping.adjPval;'])
        
    %prepare inputs for the rep mets function
    genes           = DEmapping.modelGenes;
    genePvalues     = DEmapping.adjPval;
    geneFoldChanges = DEmapping.log2FC;
    
    %Run reporter metabolites analysis!
    outputFile      = ['../../results/reporter_metabolites/repMets_glu_vs_' str '.txt'];
    repMets=get_repMets(model,genes,genePvalues,toExclude,true,outputFile,geneFoldChanges);
end    