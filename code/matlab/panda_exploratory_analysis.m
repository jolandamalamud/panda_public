%% Exploratory statistical analyses investigating the gonogo task parameter in the PANDA RCT
% 1 = placebo, 2 = sertraline
addpath(genpath(pwd));
filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/';
load([filepath 'gng_preprocessed/rctdata.mat']);
load([filepath 'gng_preprocessed/gngdata.mat']);
extract_parameters;

Nsj = length(rctdata.subid); % number of patients
Np = size(gng.params, 1); % number of cognitive parameters
T = 3; % number of sessions

% table of stratification variables
stratification_tbl = array2table(rctdata.strat_variables, ...
    'VariableNames', rctdata.strat_variables_names);
% table of baseline variables (age, sex)
baseline_variables_tbl = array2table(rctdata.baseline_variables(:,[2,9]), ...
    'VariableNames', rctdata.baseline_variables_names([2,9]));

% load raw gonogo task data
D = importdata([filepath 'panda_data/panda_datafile_date.mat']);
% D = D(vertcat(D.sess) == '1' | vertcat(D.sess) == '2');
% basic characteristics of gonogo task data (e.g. go-probability, performance)
bc = gng_basic_characteristics(D);

group = rctdata.group_allocation==2;
rctdata.sessions = [0 2 6 12]; %[1:4];

% initialize plotting class
vap = VisualAnalysisPANDA;
% score over time split in groups
% vap.score_over_time(rctdata.logtotgad, group)

% match gng and rct data based on subject IDs
[id_matching] = match_id([1:length(gng.subid)], gng.sess, rctdata.subid, gng.subid);
[gng_pcorr] = match_id(bc.pcorr, gng.sess, rctdata.subid, gng.subid);

%% Analyse basic task characteristics (all task runs included)
% different task measures: performance for each condition, overall
% performance, g2w minus ng2w and ng2a minus g2a, and reaction times
task_measure_labels = [vap.task_condition_labels, 'overall', ...
        'g2w-ng2w', 'ng2a-g2a', 'rt'];
for k = 1:3
    task_measures_sessions{k} = array2table([sq(gng_pcorr(:,:,k))', ...
    sq(nanmean(gng_pcorr(:,:,k),1))', G2WvsNG2W(:,k), ...
    NG2AvsG2A(:,k), rts(:,k)], 'VariableNames', task_measure_labels);
end

% relation between performance and psychiatric score at baseline
score_name = {'gad', 'phq'};
for k = 1:2
    task_measures = task_measures_sessions{1};
    eval(['score = rctdata.logtot' score_name{k} '(:,1);'])
    for i = 1:width(task_measures)
        [c(i),p(i)] = corr(table2array(task_measures(:,i)), ...
            score, 'rows', 'complete');
    end
    disp(['correlation between performance and ' score_name{k} 'score at baseline:'])
    disp(array2table([c',p'],'VariableNames', {'corrcoef', 'pval'}, ...
        'RowNames', task_measure_labels))
    disp('---------------------------------------------------------------')
    clear c p
end

[c, p] = gng_pearson_corr(rctdata.logtotgad, bc.pcorr, gng, rctdata, opt);

% difference of performance between groups
doprint = true;
opt.match_data = true; opt.exlude_data = true;

task_characteristics = [bc.pcorr; nanmean(bc.pcorr,1); bc.pcorr_pav; ...
bc.mrt; sq(nanmean(bc.switch_incondition,1)); sq(nanmean(bc.switch_incondition_after_worse,1))];
[results.performance.mean, t, p] = ...
    gng_simple_ttest_group_difference(task_characteristics, gng, rctdata, opt);

results.performance.stats = gng_results_table({t(:,2:3),p(:,2:3)}, ...
    {{'tstats','pval'}, {'FU1','FU2'}}, [vap.task_condition_labels, ...
     {'overall', 'G2WvsNG2W','NG2AvsG2A', 'rt', 'switching', ...
    'switching after worse outcome'}], doprint, ...
    'group difference in basic task characteristics');

%% Analyse modelling characteristics
options.disp = false;

% match gng and rct data based on subject IDs
[id_matching] = match_id([1:length(gng.subid)], gng.sess, rctdata.subid, gng.subid);
[gng_parameter] = match_id(gng.params, gng.sess, rctdata.subid, gng.subid);
model_fit = match_id(gng.modelfit, gng.sess, rctdata.subid, gng.subid);
[missing_factor] = match_id(gng.ex, gng.sess, rctdata.subid, gng.subid);
missing_factor(isnan(missing_factor)) = 2;
% 0=not missing; 1=missing due to uinformative task data; 2=missing task data

% exclude uninformative task runs
[gng_parameter_included] = exclude_uninformative_data(gng_parameter, missing_factor > 0);
% [gng_pcorr_included] = exclude_uninformative_data(gng_pcorr, missing_factor > 0);
% [raw_behaviour_included] = exclude_uninformative_data(raw_behaviour, missing_factor > 0);
% [rts_included] = exclude_uninformative_data(rts, missing_factor > 0);

% parameter over time
vap.data_over_time(gng_parameter_included);

%% any parameters (other than aversive Pav) related group allocation?
options.group = rctdata.group_allocation == 2;
options.time = [0, 2, 6];
options.confounders = stratification_tbl;
options.random_slopes = false;
options.group_time_interaction = false;
group_over_time = [zeros(Nsj,1), repmat(options.group, 1, 3-1)];

for p = 1:Np
    cog_parameter = sq(gng_parameter_included(p,:,:));
    % over all 3 sessions [0,2,6]
    options.labels = {vap.parameter_label{p}, 'group'}; options.T = 3;
    lmeT3 = run_lme(options, cog_parameter, group_over_time, ...
        options.labels, options.time, stratification_tbl);
    % over all 3 sessions [0,2,6] with group-time interaction
	options.group_time_interaction = true;
    lmeT3int = run_lme(options, cog_parameter, group_over_time, ...
        options.labels, options.time, stratification_tbl);
    % over all 2 sessions [0,2] and [0,6]
    options.T = 2; options.group_time_interaction = false;
    results = lme_fu_separate(cog_parameter, [], options);
    rr_group(:,:,p) = [results{1}.Coefficients.Estimate(2), ...
        results{1}.Coefficients.pValue(2); ....
        results{2}.Coefficients.Estimate(2), ...
        results{2}.Coefficients.pValue(2); lmeT3.Coefficients.Estimate(2), ...
        lmeT3.Coefficients.pValue(2); lmeT3int.Coefficients.Estimate(2), ...
        lmeT3int.Coefficients.pValue(2)];
end
    
stats.parameter_group = rr_group;
vap.data_over_time_groupsplit2(gng_parameter_included, group, ...
    {'placebo','sertraline'}, rr_group);    
disp('--------------------------------------------------------')

%% any parameters related to log-transformed total gad score?
T = 3;
group = rctdata.group_allocation == 2; % group allocation (1 = sertraline, 0 = placebo)
group_over_time = [zeros(Nsj,1), repmat(group,1,T-1)];

stratification_tbl = array2table(rctdata.strat_variables, ...
    'VariableNames', rctdata.strat_variables_names);

baseline_variables_tbl = array2table(rctdata.baseline_variables(:,[2,3,9]), ...
    'VariableNames', rctdata.baseline_variables_names([2,3,9]));

% run LME
options.T = T;
options.random_slopes = true;
options.group_time_interaction = false;
rr_gad = zeros(2,8);
for p = 1:8
    independent_variables = cat(3, ...
        sq(gng_parameter_included(p,:,1:3)), group_over_time);
    lme = run_lme(options, rctdata.logtotgad(:,1:3), ....
        independent_variables, {'anxiety', vap.parameter_label{p}, ...
        'group'}, [0,2,6], [stratification_tbl, baseline_variables_tbl]);
    rr_gad(:,p) = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
end

tab = array2table(rr_gad, 'VariableNames', vap.parameter_label, ...
    'RowNames', {'effect','pval'});
stats.parameter_gad = tab;
disp(tab); clear tab;
disp('higher anxiety is related to reduced loss sensitivity and higher loss lr!')

%% PCA on GAD-7 questionnaire items
firstpc_gad = dim_reduction_3d(rctdata.gaditems, 1);

% correlation log-transformed total gad score and first pc loadings
firstpc_gad_long = sq(firstpc_gad)';
pc_total_corr = corr(firstpc_gad_long(:), rctdata.logtotgad(:), ...
    'rows', 'complete');
disp(['correlation between total gad score and first pc loadings: ' ...
    num2str(round(pc_total_corr,2))]);

%% PCA on PHQ-9 questionnaire items
firstpc_phq = dim_reduction_3d(rctdata.phqitems, 1);

% correlation log-transformed total phq score and first pc loadings
firstpc_phq_long = sq(firstpc_phq)';
pc_total_corr = corr(firstpc_phq_long(:), rctdata.logtotphq(:), ...
    'rows', 'complete');
disp(['correlation between total phq score and first pc loadings: ' ...
    num2str(round(pc_total_corr,2))]);

%% What is going on with loss parameters and anxiety?
% split into high and low anxiety
anxiety_split = nan(Nsj,T);
for t = 1:3
    tertiles = quantile(rctdata.logtotgad(:,t),[0.325 0.675]);
    anxiety_split(rctdata.logtotgad(:,t) <= tertiles(1), t) = 0;
    anxiety_split(rctdata.logtotgad(:,t) >= tertiles(2), t) = 1;
end

% plot raw behaviour split in low and high anxiety group
vap.raw_behaviour_groupsplit2(raw_behaviour_included, anxiety_split, ...
    missing_factor, {'low anxiety','high anxiety'});
% plot loss parameter split in low & high anxiety and treatment group
vap.data_over_time_groupsplit4(gng_parameter_included, group, ...
    anxiety_split, {'low anxiety','low anxiety 5HT','high anxiety', ...
    'high anxiety 5HT'});

%% PCA on loss LR and loss sensitivity
% correlation between parameters
[paramter_corr, pval] = corr(reshape(gng_parameter_included, Np, ...
    Nsj*T)', 'rows', 'complete');
stats.parameter_correlation.corr = tril(pval) + triu(paramter_corr);
stats.parameter_correlation.alpha = 0.05 / (Np * (Np-1) / 2);

disp(['loss lR and loss sensitivity are highly negatively correlated: ' ...
    num2str(stats.parameter_correlation.corr(2,4)), ' , p = ' ...
    num2str(stats.parameter_correlation.corr(4,2))])

[pcloadings, pccomponents] = ...
    dim_reduction_3d(gng_parameter_included([2,4],:,:), 1);
disp(['first principle component: ' num2str(pccomponents(:,1)')])

T = 3;
group_over_time = [zeros(Nsj,1), repmat(group,1,T-1)];
independent_variables = cat(3,sq(pcloadings(1,:,1:T)), ...
    group_over_time);

% run LME
lme = run_lme(options, rctdata.logtotgad(:,1:3), independent_variables, ...
    {'anxiety', 'pcloadings_loss', 'group'}, [0,2,6], ...
    [stratification_tbl, baseline_variables_tbl]);

% save results
stats.parameter_pc_loss.lme = lme;
stats.parameter_pc_loss.results = ...
    [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)];

%% PCA on all parameter
disp(array2table(stats.parameter_correlation.corr, 'VariableNames', ...
    vap.parameter_label));
disp(['multiple comparison: ' num2str(stats.parameter_correlation.alpha)]);

% PCA on parameter - participants loadings on the first 3 PCs
[pcloadings, pcomponents] = dim_reduction_3d(gng_parameter_included, 3);
plot_title = {'1st PC: loss lr VS loss sensitivity, aversive Pav and noise', ...
    '2nd PC: reward lr and go bias VS noise', ...
    '3rd PC: reward lr and Pav biases VS reward sensitivity, go bias and noise'};
vap.plot_pca(pcomponents(:,1:3), plot_title);

%% Does change predict treatment outcome?
[par_change] = slope(gng_parameter_included(:,:,1:3));

rr_gad_par_slope = nan(8,2);
gad0 = rctdata.logtotgad(:,1);
gad4 = rctdata.logtotgad(:,3);
for k = 1:8
    par_slope = sq(par_change(1,k,:));
    lm = fitlm(table(gad4,gad0, par_slope, group), 'gad4 ~ gad0 + par_slope + group');
    rr_gad_par_slope(k,:) = [lm.Coefficients.Estimate(3), lm.Coefficients.pValue(3)];
end
stats.slopes_parameter_gad4 = rr_gad_par_slope;
disp(rr_gad_par_slope);
disp('beta of change in go bias predicting anxiety at next time point is lower in sertraline group');

p = nan(2,1); c = nan(2,1);
for i = 1:2; gg= i-1; [c(i),p(i)] = corr(gad_slope(group==gg), ...
    sq(par_change(1,8,group==gg)), 'rows', 'complete');
end
disp(['corr in placebo group: ' num2str(c(1)) ' , p = ' num2str(p(1)) ...
    ', in sertraline group: ' num2str(c(2)) ' , p = ' num2str(p(2))]);

figure(); 
subplot(121); plot(sq(par_change(1,8,group==0)), gad4(group==0), 'ro'); lsline
title('placebo group');
xlabel('change of go bias over sessions'); ylabel('anxiety at 12 weeks');
subplot(122); plot(sq(par_change(1,8,group==1)), gad4(group==1), 'bo'); lsline
xlabel('change of go bias over sessions'); ylabel('anxiety at 12 weeks');
title('sertraline group');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% Task performance and sertraline, resp. psychiatric scores

rr_performance_group = zeros(2,4);
group = rctdata.group_allocation == 2;
for i = 1:size(gng_pcorr_included,1)
    performance = sq(gng_pcorr_included(i,:,1:2));
    lme = run_lme(performance, [zeros(Nsj,1), group], {'performance','group'}, false, [0, 2 ,6]);
    rr_performance_group(:,i) = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
end
disp('probability correct related to group?')
tab = array2table(rr_performance_group, 'VariableNames', vap.task_condition_labels, ...
    'RowNames', {'effect','pval'});
stats.performance.group = tab;
disp(tab); clear tab;
disp(['multiple comparison p < ' num2str(0.05 / 4)])

rr_performance_gad = zeros(2,4);
for i = 1:size(gng_pcorr_included,1)
    performance = sq(gng_pcorr_included(i,:,:));
    lme = run_lme(rctdata.logtotgad(:,1:3), performance, [zeros(Nsj,1), group, group]);
    rr_performance_gad(:,i) = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
end
disp('probability correct related to anxiety?')
tab = array2table(rr_performance_gad, 'VariableNames', vap.task_condition_labels, ...
    'RowNames', {'effect','pval'});
stats.performance.gad = tab;
disp(tab); clear tab;
disp(['multiple comparison p < ' num2str(0.05 / 4)])

%% Reaction Times
% reaction time related to performance?
disp(corr(nanmean(gng_pcorr_included,3)', ...
    sq(nanmean(rts_included,2)), 'rows', 'complete'));
disp('the faster the worse in go condition, vice versa in nogo condition');

% reaction time related to cognitive parameters?
for k = 1:8
    disp(corr(nanmean(gng_parameter_included(k,:,:),3)', ...
        sq(nanmean(rts2_included,2)),'rows', 'complete')); 
end

% reaction time related to group?
vap.data_over_time_groupsplit2(sq(nanmean(rts_included,1)), group, ...
    {'placebo','5HT'})

%% cognitive parameters related to baseline measures
for i = 1:length(rctdata.strat_variables_names)
    for k = 1:3
    [c,p]=corr(rctdata.strat_variables(:,i), sq(gng_parameter(:,:,k))', ...
        'rows', 'complete'); 
    tt = [c;p]; disp(rctdata.strat_variables_names{i}); disp(array2table(round(tt(:),2)'));
    end
end

%% parameter related to single log-transformed psychiatric questionnaire items

% gad
rr_gad = zeros(2,8,7);
for p = 1:8
    for i = 1:7
        singleItem_score = sq(rctdata.gaditems(i,1:3,:))';
        parameter = sq(gng_parameter_included(p,:,:));
        lme = gng_lme_control_baseline(singleItem_score, parameter, rctdata);
        rr_gad(:,p,i) = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
    end
end

rownames = [];for i = 1:7; rownames = [rownames,{['effect gad-item ' ...
        num2str(i)], ['pval gad-item ' num2str(i)]}]; end
tab = array2table(reshape(rr_gad, 14,8), 'VariableNames', ...
    vap.parameter_label, 'RowNames', rownames);
stats.subitems.gad = tab; 
disp(tab); clear tab;
disp(['multiple comparison p < ' num2str(0.05 / (7*8))])

% phq
rr_phq = zeros(2,8,9);
for p = 1:8
    for i = 1:9
        singleItem_score = sq(rctdata.phqitems(i,1:3,:))';
        parameter = sq(gng_parameter_included(p,:,:));
        lme = gng_lme_control_baseline(singleItem_score, parameter, rctdata);
        rr_phq(:,p,i) = [lme.Coefficients.Estimate(2),lme.Coefficients.pValue(2)]; 
    end
end

rownames = [];for i = 1:9; rownames = [rownames,{['effect phq-item ' ...
        num2str(i)], ['pval phq-item ' num2str(i)]}]; end
tab = array2table(reshape(rr_phq, 18, 8), 'VariableNames', ...
    vap.parameter_label, 'RowNames', rownames);
stats.subitems.phq = tab; 
disp(tab); clear tab;
disp(['multiple comparison p < ' num2str(0.05 / (9*8))])

%% save statistical results
goalpath = '~/phd/projects/gng_panda_antler/gng_panda/results/';
save([goalpath, 'gng_exploratory_results'], 'stats'); clear all;