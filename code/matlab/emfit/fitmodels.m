addpath(genpath(pwd));
filepath = '../../';
% loaddata
D = importdata([filepath 'data/panda_data/panda_datafile_date.mat']);
resultsDir = [filepath 'results/model_fits/'];

mkdir(resultsDir)
data = D;
bc = gng_basic_characteristics(data);

% Define data to fit
% for t = 1:3
%     resultsDir = [filepath 'results/model_fits/sess' num2str(t)];

% data = data(nanmean(bc.pcorr,1) >= 0.5);
% combine participants sessions
% data = combine_sessions(D);
% bc = gng_basic_characteristics(data);

rng('default');


% FITTING USING EM
modelsToFit = [1,3,4,5,6,8];
models=modelList;
models = models(modelsToFit);
batchRunEMfit('mAffectiveGoNogo', data, resultsDir, 'modelstofit', ...
    modelsToFit, 'checkgradients', 1);
% end