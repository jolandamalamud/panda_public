filepath = '~/phd/projects/antler/';
resultsDir = [filepath 'results/recovery'];

%% Run parameter recovery for the different models
modelfiles = dir([filepath '/ll*.mat']);

for f = 1:length(modelfiles)
    model.name = modelfiles(f).name(1:end-4);

    R = load([filepath '/' model.name '.mat']);
    model.npar = size(R.E,1);

    [Data, fit] = gng_generate_data_and_emfit(model, R, 1000, 96);
     
    % save fit
    save([resultsDir '/' model.name '_recovery.mat'], 'model', 'Data', 'fit');
end

% mse and correlation between true and estimated parameters
for f = 1:length(modelfiles)
    model.name = modelfiles(f).name(1:end-4);
    R = load([resultsDir '/' model.name '_recovery.mat']);
    R = gng_parameter_mse_corr(R);
    rec.(model.name) = R;
    plot_data.partrue(1:R.model.npar,:,f) = R.trueparams;
    plot_data.parest(1:R.model.npar,:,f) = R.fit.parest;
    plot_data.pcc(f,1:R.model.npar) = R.pcc;
end

modellist;
gng_plot_parameter_recovery(plot_data, model_characteristics, 1);

%% run parameter recovery multiple times for the winning model
ms = 100; 
model.name = 'll2b2a2epxb';
R = load([filepath '/data/modelling_results/' model.name '.mat']);
modelsToFit = [8];
models = modelList;
models = models(modelsToFit);

for f = 1:ms
    [Data, modelfit] = gng_generate_data_and_emfit(models, R, 200, 96);
    D{f} = Data;
    fit{f} = modelfit; 
end
tv = datestr(now, 'yyyymmddHHMM');
save([resultsDir '/' model.name '_recovery_' tv '.mat'], 'model', 'D', 'fit', '-v7.3');
    
% mse and correlation between true and estimated parameters
for f = 1:ms
    R.Data = D{f}; R.fit = fit{f};
    R.model = models;
    R = gng_parameter_mse_corr(R);
    plot_pcc(:,f) = R.pcc;
end

gng_plot_parameter_recovery(plot_pcc, model, 0);

%%
%% run parameter recovery multiple times for the winning model
ms = 10; 
model.name = 'll2b2a2epxb';
R = load([filepath '/results/model_fits/' model.name '.mat']);
modelsToFit = [8];
models = modelList;
models = models(modelsToFit);

% used only informative data for recovery
R2 = load([filepath '/results/model_fits/llb.mat']);
ex = R.bf.iL - R2.bf.iL < 3;
R.E = R.E(:,~ex);

parfor f = 1:ms
    [Data, modelfit] = gng_generate_data_empirical_and_fit(models, R, sum(ex==0), 96);
    D{f} = Data;
    fit{f} = modelfit; 
end
tv = datestr(now, 'yyyymmddHHMM');
save([resultsDir '/' model.name '_recovery_' tv '.mat'], 'model', 'D', 'fit');
    
% % mse and correlation between true and estimated parameters
% for f = 1:ms
%     R.Data = D{f}; R.fit = fit{f};
%     R.model = models;
%     R = gng_parameter_mse_corr(R);
%     plot_pcc(:,f) = R.pcc;
% end
% 
% gng_plot_parameter_recovery(plot_pcc, model, 0);