% independent session modelling VS parameter change modelling
addpath(genpath(pwd));
filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/';
load([filepath 'gng_preprocessed/gngdata.mat']);
load([filepath '/modelling_results/ll2b2a2epxbdbxbl.mat'])
gng_parameter = match_id(gng.params, gng.sess, unique(gng.subid), gng.subid);

vap = VisualAnalysisPANDA;

% correlation with baseline parameter estimates
for k = 1:8
    [c(k),p(k)] = corr(E(k, :)', gng_parameter(k,:,1)', 'rows', 'complete'); 
end

% correlation with session 1 and 2 and change parameter estimates
ii = [8, 7, 2]; iii = [9, 10, 11];
for i = 1:3
    for k = 1:2
        [c(8+k*i+(k-2)*(1-i)),p(8+k*i+(k-2)*(1-i))] = corr(E(ii(i), :)' + E(iii(k), :)', ...
            gng_parameter(ii(i),:,k)', 'rows', 'complete'); 
    end
end

array2table(round([c',p'],2), 'rownames', [vap.parameter_label, ...
    'delta_bias1', 'delta_bias2', 'delta_noise1', 'delta_noise2', ...
    'delta_sensL1', 'delta_sensL2'], 'variablenames', ...
    {'correlation', 'pval'})