function [Data, fit] = gng_generate_data_and_emfit(model, R, Nsj, Ntrials)
    
    % take the estimated hyperparameters to generate data
    setparams.sigma = R.stats.groupvar;
    setparams.mu = mean(R.stats.subjectmeans,2);

    % draw parameters from hyperprior
    params = mvnrnd(setparams.mu, setparams.sigma, Nsj);
    params = R.E';
    % generate data
    Data = generateExampleDataset(Nsj, [], model, params', Ntrials);

    %% fit simulated data using EM
    fit.parest = emfit(model.name, Data, model.npar);
%     [fit.parest,fit.V,fit.alpha,fit.stats,fit.bf] = emfit(model.name, Data, model.npar);

    
    % orthogonalize parameters 
    % Q = orth(params);
    % Data_orth=generateExampleDataset(1000,resultsDir, model, Q', 96);
    % [fitorth.parest] = emfit(model.name,Data_orth,model.npar);\
    
end