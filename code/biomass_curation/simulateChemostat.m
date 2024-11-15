function mod_data = simulateChemostat(model,exp_data,GAM,POratio,NGAM)
%Modify GAM withouth changing the protein content:
parameters.exch_names{1} = 'growth';
parameters.exch_names{2} = 'D-glucose exchange';
parameters.exch_names{3} = 'oxygen exchange';
parameters.exch_names{4} = 'carbon dioxide exchange';
parameters.bioRxn = 'r_4041';
%Relevant positions:
exch_names  = parameters.exch_names;
pos(1)      = find(strcmp(model.rxnNames,exch_names{1}));
pos(2)      = find(strcmp(model.rxnNames,exch_names{2}));
pos(3)      = find(strcmp(model.rxnNames,exch_names{3}));
pos(4)      = find(strcmp(model.rxnNames,exch_names{4}));
if nargin>4
    model = changeNGAM(model,NGAM);
elseif nargin>3
    model = changePOratio(model,POratio);
elseif nargin>2 
    model = changeGAM(model,GAM);
end
%Simulate chemostats:
mod_data = zeros(size(exp_data));
for i = 1:length(exp_data(:,1))
    %Fix biomass
    model = setParam(model,'lb',model.rxns(pos(1)),exp_data(i,1));
    %set an arbitrarily high glucose uptake rate
    model = setParam(model,'ub',model.rxns(pos(2)),0);
    model = setParam(model,'lb',model.rxns(pos(2)),-10);
    %minimize glucose
    model = setParam(model,'obj',model.rxns(pos(2)),1);
    sol   = solveLP(model,1);
    %printFluxes(model,sol.x,true)
    %pause
    %Store relevant variables:
    mod_data(i,:) = sol.x(pos)';
end
end