%% Supplementary statistical analyses - Controlling for time in anlyses!!!
addpath(genpath(pwd));
filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/gng_preprocessed/';
load([filepath 'rctdata.mat']);
% load([filepath 'gngdata.mat']);
extract_parameters;

Nsj = length(rctdata.subid); % number of patients
Np = size(gng.params, 1); % number of gonogo task parameters
T = max(gng.sess); % number of task sessions

% match gng and rct data based on subject IDs
[id_matching] = match_id([1:length(gng.subid)], gng.sess, rctdata.subid, gng.subid);
[gng_parameter] = match_id(gng.params, gng.sess, rctdata.subid, gng.subid);
% gng.ex = gng.ex | nanmean(bc.pcorr,1) < 0.5;
[missing_factors] = match_id(gng.ex, gng.sess, rctdata.subid, gng.subid);
missing_factors(isnan(missing_factors)) = 2;
% 0=not missing; 1=missing due to uinformative task data; 2=missing task data

% exclude uninformative task runs
[gng_parameter_included] = exclude_uninformative_data(gng_parameter, missing_factors > 0);

%% define variables for analysis
% group allocation variable (1 = sertraline, 0 = placebo)
options.group = rctdata.group_allocation == 2;
% time variable (in weeks)
options.time = [1, 2, 3];
% table of stratification variables
stratification_tbl = array2table(rctdata.strat_variables, ...
    'VariableNames', rctdata.strat_variables_names);
% table of baseline variables (age, sex)
baseline_variables_tbl = array2table(rctdata.baseline_variables(:,[2,9]), ...
    'VariableNames', rctdata.baseline_variables_names([2,9]));
% aversive Pavlovian bias
cog_parameter = sq(gng_parameter_included(4,:,:));
% psychiatric score
psy_score = rctdata.logtotgad;
% disp all statistics
options.disp = true;
%% Hypothesis 1: Does sertraline relate to the aversive Pavlovian bias?
% Mixed-effects model correcting for time!!!
options.T = 3;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];
% run LME
options.random_slopes = false;
options.group_time_interaction = false;
lme = run_lme(options, cog_parameter, group_over_time, ...
    {'avPavBias', 'group'}, options.time, stratification_tbl)

disp('drug effect disappeared when including time (in weeks) as covariate')
disp('-------------------------------------------------------------------')

%% Analysis for two follow-up time points separately
options.T = 2;
options.confounders = stratification_tbl;
options.labels = {'avPavBias', 'group'};
% run LME
results = lme_fu_separate(cog_parameter, [], options);

%% Analysis for group-time interaction
options.T = 3;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];

% run LME
options.random_slopes = false;
options.group_time_interaction = true;
lme = run_lme(options, cog_parameter, group_over_time, ...
    {'avPavBias', 'group'}, options.time, stratification_tbl)
disp('--------------------------------------------------------')

%% Hypothesis 2: Does the aversive Pavlovian bias relate to anxiety?
% Mixed-effects model correcting for time!!!
options.T = 3;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];

% run LME
options.random_slopes = true;
options.group_time_interaction = false;
independent_variables = cat(3,cog_parameter, group_over_time);
lme = run_lme(options, psy_score(:,1:3), independent_variables, ...
    {'anxiety', 'avPavbias', 'group'}, options.time, ...
    [stratification_tbl, baseline_variables_tbl])

%% Analysis for two follow-up time points separately
options.T = 2;
options.confounders = [stratification_tbl, baseline_variables_tbl];
options.labels = {'anxiety', 'avPavbias', 'group'};
% run LME
results = lme_fu_separate(psy_score, cog_parameter, options);

%% Hypothesis 5: Does the appetitive Pavlovian bias relate to depression?
% Mixed-effects model correcting for time!!!
options.T = 3;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];
options.time = [1,2,3];
% run LME
options.random_slopes = true;
options.group_time_interaction = false;
independent_variables = cat(3, sq(gng_parameter_included(5,:,:)), group_over_time);
lme = run_lme(options, rctdata.logtotphq(:,1:3), independent_variables, ...
    {'phq9', 'appPavbias', 'group'}, options.time, ...
    [stratification_tbl, baseline_variables_tbl])

%% Analysis for two follow-up time points separately
options.T = 2;
options.confounders = [stratification_tbl, baseline_variables_tbl];
options.labels =  {'phq9', 'appPavbias', 'group'};
% run LME
results = lme_fu_separate(rctdata.logtotphq, sq(gng_parameter_included(5,:,:)), options);

%% Hypothesis 6: Does the reward sensitivity relate to anhedonia?
% Mixed-effects model correcting for time!!!
options.T = 3;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];

% run LME
options.random_slopes = true;
options.group_time_interaction = false;
independent_variables = cat(3, sq(gng_parameter_included(1,:,:)), group_over_time);
lme = run_lme(options, rctdata.loganhedonia(:,1:3), independent_variables, ...
    {'anhedonia', 'sensRew', 'group'}, options.time, ...
    [stratification_tbl, baseline_variables_tbl])