%generateCintGEM

%load SBML draft model and extract information for S. cerevisiae orthologs
%model = readCbModel('../models/candida_intermedia/Candida_intermedia.xml');
%git clone --depth 1 https://github.com/IVANDOMENZAIN/diverseYeasts_metabolism.git
%model = importModel('diverseYeasts_metabolism/models/candida_intermedia/Candida_intermedia.mat');
sce_proteins = model.proteins;
model = ravenCobraWrapper(model);
model.proteins = sce_proteins;
model = curateLeloirPathway(model);

%integrate galactose oxido-reductive pathway
cd oxRed_path_addition/
model = oxidopathwayaddition(model);
model = mapGeneIDs(model);
model = lxr4_exploration(model);
cd ..
model = correct_gene_annotation(model);
cd biomass_curation/
model = adjust_biomass_comp(model);
cd ..
save('../model/sintGEM.mat','model')
%generate version-controllable files
formulas = constructEquations(model);
rxns = model.rxns;
rxnNames = model.rxnNames;
grRules = model.grRules;
modelTable = table(rxns,rxnNames,formulas, grRules);
writetable(modelTable,'../model/sintGEM.txt','WriteVariableNames',true,'Delimiter','\t','QuoteStrings',false);
%add version control
genes = model.genes;
shortnames = model.geneShortNames;
orthologues = model.orthologues;
proteins = model.proteins;
gene_table = table(genes,shortnames,orthologues,proteins);
writetable(gene_table,'../model/gene_table_sintGEM.txt','Delimiter','\t','QuoteStrings',false);
%write SBML
exportModel(model,'../model/sintGEM.xml');