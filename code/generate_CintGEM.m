%generateCintGEM

%load SBML draft model and extract information for S. cerevisiae orthologs
%model = readCbModel('../models/candida_intermedia/Candida_intermedia.xml');

load('../models/candida_intermedia/Candida_intermedia.mat')
sce_proteins = model.proteins;
model = ravenCobraWrapper(model);
model.proteins = sce_proteins;
model = curateLeloirPathway(model);
load('../models/candida_intermedia/cint_leloir.mat')

%integrate galactose oxido-reductive pathway
cd oxRed_path_addition/
model = oxidopathwayaddition(model);



mapGeneIDs
lxr4_exploration
cd ..
correct_gene_annotation
cd biomass_curation/
model = adjust_biomass_comp;
cd ..
save('../models/candida_intermedia/cintGEM.mat','model')
%generate version-controllable files
formulas = constructEquations(model);
rxns = model.rxns;
rxnNames = model.rxnNames;
grRules = model.grRules;
modelTable = table(rxns,rxnNames,formulas, grRules);
writetable(modelTable,'../models/candida_intermedia/cintGEM.txt','WriteVariableNames',true,'Delimiter','\t','QuoteStrings',false);
%add version control
genes = model.genes;
shortnames = model.geneShortNames;
orthologues = model.orthologues;
proteins = model.proteins;
gene_table = table(genes,shortnames,orthologues,proteins);
writetable(gene_table,'../models/candida_intermedia/gene_table_cintGEM.txt','Delimiter','\t','QuoteStrings',false);
%remove temporary .mat files
delete '..'/models/candida_intermedia/cintGEM_curated2.mat
delete '..'/models/candida_intermedia/cintGEM_curated.mat
delete '..'/models/candida_intermedia/cint_GEM_curated.mat
delete '..'/models/candida_intermedia/cintGEM_gene_curated.mat
delete '..'/models/candida_intermedia/cintGEM_oxido.mat
delete '..'/models/candida_intermedia/cintGEM_oxido_orthologs.mat
delete '..'/models/candida_intermedia/cintGEM_oxido_orthologs_curated.mat