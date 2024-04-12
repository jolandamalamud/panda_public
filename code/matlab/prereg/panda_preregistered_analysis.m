%% Main statistical analyses investigating the gonogo task parameter in the PANDA RCT
% 1 = placebo, 2 = sertraline
addpath(genpath(pwd));
filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/gng_preprocessed/';
load([filepath 'rctdata.mat']);
load([filepath 'gngdata.mat']);

Nsj = length(rctdata.subid); % number of patients
Np = size(gng.params, 1); % number of gonogo task parameters
T = max(gng.sess); % number of task sessions

% match gng and rct data based on subject IDs
[id_matching] = match_id([1:length(gng.subid)], gng.sess, rctdata.subid, gng.subid);
[gng_parameter] = match_id(gng.params, gng.sess, rctdata.subid, gng.subid);
[missing_factors] = match_id(gng.ex, gng.sess, rctdata.subid, gng.subid);
missing_factors(isnan(missing_factors)) = 2;
% 0=not missing; 1=missing due to uinformative task data; 2=missing task data

% exclude uninformative task runs
[gng_parameter_included] = exclude_uninformative_data(gng_parameter, missing_factors > 0);

%% Missing data

% investigating missingness due to uninformative task behaviour 
missing_factor = missing_factors;
missing_factor(missing_factor==2)=nan;

% does the amount of missing data differ between group?
[h, pval] = ttest2(missing_factor(rctdata.group_allocation==1,:), ...
    missing_factor(rctdata.group_allocation==2,:), 0.05/T);
stats.missing.group = table(h', pval', ...
    nanmean(missing_factor(rctdata.group_allocation==1,:),1)', ...
    nanmean(missing_factor(rctdata.group_allocation==2,:),1)');
stats.missing.group.Properties.VariableNames = ...
    {'sig', 'pval', 'placebo', 'sertraline'};
stats.missing.group.Properties.RowNames = {'T1', 'T2', 'T3'};
disp(stats.missing.group); clear m1 m2 pval

% are baseline variables associated with missing data?
tab = ttest_missingness([rctdata.strat_variables, ...
    rctdata.baseline_variables], missing_factor);
tab.Properties.RowNames = [rctdata.strat_variables_names, ...
    rctdata.baseline_variables_names];
stats.missing.baseline = tab; 
disp(stats.missing.baseline); clear tab

% are gad7 and phq9 associated with missing data?
tab = ttest_missingness(cat(3,rctdata.logtotgad(:,1:3), ...
    rctdata.logtotphq(:,1:3)), missing_factor);
tab.Properties.RowNames = {'GAD7', 'PHQ9'};
stats.missing.phqgad = tab; 
disp(stats.missing.phqgad); clear tab

%% Hypothesis 1: Does sertraline relate to the aversive Pavlovian bias?
% Mixed-effects model
group = rctdata.group_allocation == 2; % group allocation (1 = sertraline, 0 = placebo) 
neg_pav_all = sq(gng_parameter_included(6,:,:)); % aversive Pavlovian parameter
cistot = rctdata.strat_variables(:,1); % stratification variables
depdur = rctdata.strat_variables(:,2); % stratification variables
site = rctdata.strat_variables(:,3); % stratification variables

% transform variables for LME
subject = repmat(1:Nsj, T, 1); subject = subject(:); 
neg_pav_lme = reshape(neg_pav_all', Nsj*T, 1);
group_lme = [zeros(Nsj,1), group, group];
group_lme = reshape(group_lme', Nsj*T, 1);
cistot_lme = repmat(cistot',T,1); cistot_lme = cistot_lme(:); 
depdur_lme = repmat(depdur',T,1); depdur_lme = depdur_lme(:); 
site_lme = repmat(site',T,1); site_lme = site_lme(:);
% run LME
lme = fitlme(table(neg_pav_lme, group_lme, cistot_lme, depdur_lme, ...
    site_lme, subject), ['neg_pav_lme ~ group_lme + cistot_lme + ' ... 
    'depdur_lme + site_lme + (1|subject)']);

% save results
stats.H1.lme = lme;
stats.H1.results = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
disp(lme);

%% Hypothesis 2: Does the aversive Pavlovian bias relate to anxiety?
% Mixed-effects model
gad_all = rctdata.logtotgad(:,1:3); % GAD7 total score
age = rctdata.baseline_variables(:,2); % confounders (baseline variables) 
education = rctdata.baseline_variables(:,4); % confounders (baseline variables) 
AD_past = rctdata.baseline_variables(:,5); % confounders (baseline variables)

% transform variables for LME
subject = repmat(1:Nsj,T,1); subject = subject(:);
gad_lme = reshape(gad_all', Nsj*T ,1);
age_lme = repmat(age',T,1); age_lme = age_lme(:);
education_lme = repmat(education',T,1); education_lme = education_lme(:); 
AD_past_lme = repmat(AD_past',T,1); AD_past_lme = AD_past_lme(:);

% run LME
lme = fitlme(table(gad_lme, neg_pav_lme, group_lme, cistot_lme, ...
depdur_lme, site_lme, age_lme, education_lme, AD_past_lme, subject), ...
['gad_lme ~ neg_pav_lme + group_lme + cistot_lme + depdur_lme + ' ... 
'site_lme + age_lme + education_lme + AD_past_lme + ' ... 
'(neg_pav_lme-1|subject) + (1|subject)']);

% save results
stats.H2.lme = lme;
stats.H2.results = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
disp(lme);

%% Hypothesis 3: Does the aversive Pavlovian bias mediate the effect of sertraline on anxiety?
% Multilevel mediation

% only if H1 and H2 are significant
if stats.H1.results(2) <= 0.05 && stats.H2.results(2) <= 0.05

% set up bootstrapping
opt = statset('UseParallel',true);
n_sampling = 10000;

ff1 = @(neg_pav_lme, group_lme, subject)retrieve_estimate_lme(...
    table(neg_pav_lme, group_lme, subject), ...
    'neg_pav_lme ~ group_lme + (1|subject)');
bootstat1 = bootstrp(n_sampling, ff1, neg_pav_lme, group_lme, subject, 'Options', opt);
stats.H3.a = mean(bootstat1(:,2));


ff2 = @(gad_lme, neg_pav_lme, group_lme, subject)retrieve_estimate_lme(...
    table(gad_lme, neg_pav_lme, group_lme, subject), ...
    'gad_lme ~ group_lme + neg_pav_lme + group_lme*neg_pav_lme + (neg_pav_lme-1|subject) + (1|subject)');
bootstat2 = bootstrp(n_sampling, ff2, gad_lme, neg_pav_lme, group_lme, subject, 'Options', opt);
stats.H3.b = mean(bootstat2(:,3));
stats.H3.th = mean(bootstat2(:,4));

% pure natural indirect effect -> in control group
stats.H3.pnie = mean(bootstat1(:,2) .* bootstat2(:,3));
stats.H3.ci_ab_pnie = calculate_CI(bootstat1(:,2) .* bootstat2(:,3), ...
    95, true); % 95% confidence interval
stats.H3.sig_pnie = ~(0 >= results.ci_ab_pnie(1) & 0 <= results.ci_ab_pnie(2));

% total natural indirect effect -> in treatment group
stats.H3.tnie = mean(bootstat1(:,2) .* bootstat2(:,3) + ...
    bootstat1(:,2) .* bootstat2(:,4));
stats.H3.ci_ab_tnie = calculate_CI(bootstat1(:,2) .* bootstat2(:,3) ...
    + bootstat1(:,2) .* bootstat2(:,4), 95, true); % 95% confidence interval
stats.H3.sig_tnie = ~(0 >= results.ci_ab_tnie(1) & 0 <= results.ci_ab_tnie(2));

% check direction of effect - time-lagged mixed-effect regression
gad12 = rctdata.gad7(:,4); % GAD7 total score at FU3 (12 weeks)
gad12_lme = repmat(gad12',T,1); gad12_lme = gad12_lme(:);
gad0 = rctdata.gad7(:,1); % baseline GAD7
gad0_lme = repmat(gad0',T,1); gad0_lme = gad0_lme(:);
lme = fitlme(table(gad12_lme, neg_pav_lme, gad0_lme, group_lme, subject),...
    'gad12_lme ~ neg_pav_lme + gad0_lme + group_lme + (neg_pav_lme-1|subject) + (1|subject)');

stats.H3.lme = lme;
stats.H3.results = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)];
disp(lme);
end

%% Hypothesis 4: Does the baseline aversive Pavlovian bias predict future anxiety score in the treatment group?
% time-lagged inear regression

neg_pav0 = sq(gng_parameter_included(6,:,1))'; % baseline aversive Pavlovian bias
gad0 = rctdata.logtotgad(:,1);
gad12 = rctdata.logtotgad(:,4);

lm = fitlm(table(gad12, neg_pav0, group, gad0), ...
    'gad12 ~ neg_pav0 + group + gad0 + group*neg_pav0');

stats.H4.lm = lm;
stats.H4.results = [lm.Coefficients.Estimate(5),lm.Coefficients.pValue(5)];
disp(lm);

%% Hypothesis 5: Does the appetitive Pavlovian bias relate to depression?
% Mixed-effects model

phq_all = (rctdata.logtotphq(:,1:3)); % PHQ9 total score
pos_pav_all = sq(gng_parameter_included(5,:,:)); % appetitive Pavlovian bias

phq_lme = reshape(phq_all', Nsj*T ,1);
pos_pav_lme = reshape(pos_pav_all', Nsj*T ,1);

lme = fitlme(table(phq_lme, pos_pav_lme, group_lme, subject), ...
    'phq_lme ~ pos_pav_lme + group_lme + (pos_pav_lme-1|subject) + (1|subject)');

stats.H5.lme = lme;
stats.H5.results = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)];
disp(lme);

%% Hypothesis 6: Does the reward sensitivity relate to anhedonia?
% Mixed-effects model

anhedonia_all = (rctdata.anhedonia(:,1:3)); % anhedonia score
rew_sens_all = sq(gng_parameter_included(1,:,:)); % reward sensitivity parameter

anhedonia_lme = reshape(anhedonia_all', Nsj*T ,1);
rew_sens_lme = reshape(rew_sens_all', Nsj*T ,1);

lme = fitlme(table(anhedonia_lme, rew_sens_lme, group_lme, subject), ...
    'anhedonia_lme ~ rew_sens_lme + group_lme + (rew_sens_lme-1|subject) + (1|subject)');

stats.H6.lme = lme;
stats.H6.results = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)];
disp(lme);

%% save statistical results
goalpath = '~/phd/projects/gng_panda_antler/gng_panda/results/';
save([goalpath, 'gng_prereg_results'], 'stats'); clear all;