%% plots for panda gonogo task modelling

filepath = '/Users/jolandamalamud/phd/projects/gng_panda_antler/gng_panda/';
% load empirical data
D = importdata([filepath 'data/panda_data/panda_datafile_date.mat']);

% load surrogate data
% files = dir([filepath 'data/modelling_results/SurrogateData_*']);
% files(1).name = 'SurrogateData_llbdb.mat';
files(1).name = 'SurrogateData_llba2epxbm.mat';

for f= 1:length(files)
    surrogate_data = importdata([filepath 'data/modelling_results/sess1/' files(f).name]);
    model(f).name = files(f).name(15:end-4);

    for s = 1:size(surrogate_data,1)
        surr(s).(model(f).name) = surrogate_data(s,:);
    end
end

%% independent session modelling VS parameter change modelling
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

%% prepare data to plot
task = combine_sessions(D);
for k = 1:length(task)
    plot_data(k).a = vertcat(task(k).data.a);
    plot_data(k).r = vertcat(task(k).data.r);
    plot_data(k).s = vertcat(task(k).data.s);
end
plot_surrdata = surr;
mm = model(1);

%% Figure to check overall model performance
pp = gng_plotresultsClass; %close gcf;
[as, bs] = pp.average_go(plot_data, plot_surrdata, mm);
figure(); pp.plot_model_performance_basic(as, bs, [])

%% Figure included vs excluded subjects based on iL
ex = gng_exclude('llbdb', 'll2b2a2epxbdbxbl')';
opt.ll = '-'; opt.lw = 2; opt.ll2 = '--'; opt.lw2 = 2;
opt.cl = [1 0 0 0.5]; opt.cl2 = [1 0 0];
figure(); pp.plot_model_performance_basic(as(:,:,ex), bs(:,:,ex,1), [], opt)
opt.cl = [0 1 0 0.5]; opt.cl2 = [0 1 0];
hold on; pp.plot_model_performance_basic(as(:,:,~ex), bs(:,:,~ex,2), [], opt)
legend( {'excluded','model','included', 'model'})
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% only included subjects
mm = model(2);
[as, bs] = pp.average_go(plot_data(~ex), plot_surrdata(~ex), mm);

%% model performance in both groups
group = match_id(rctdata.group_allocation==2, ones(length(unique(rctdata.subid)),1), unique(gng.subid), rctdata.subid);
pp = gng_plotresultsClass; close gcf;
[as, bs] = pp.average_go(plot_data, plot_surrdata, mm);
opt.cl = [1 0 0 0.5]; opt.cl2 = [1 0 0];
figure(); pp.plot_model_performance_basic(as(:,:,group==0), bs(:,:,group==0), [], opt)
opt.cl = [0 0 1 0.5]; opt.cl2 = [0 0 1];
hold on; pp.plot_model_performance_basic(as(:,:,group==1), bs(:,:,group==1), [], opt)
legend( {'placebo','model placebo','sertraline', 'model sertraline'})
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% model performance in different anxiety groups
baseline_gad = match_id(rctdata.logtotgad(:,1), ones(length(unique(rctdata.subid)),1), unique(gng.subid), rctdata.subid);
baseline_gad_split = quantile(baseline_gad(~ex),2);
baseline_low_gad = baseline_gad(~ex) <= baseline_gad_split(1);
baseline_high_gad = baseline_gad(~ex) <= baseline_gad_split(2);
opt.cl = [0 1 0 0.5]; opt.cl2 = [0 1 0]; 
figure(); pp.plot_model_performance_basic(as(1:24,:,baseline_low_gad), bs(1:24,:,baseline_low_gad), [], opt)
opt.cl = [1 0 0 0.5]; opt.cl2 = [1 0 0];
hold on; pp.plot_model_performance_basic(as(1:24,:,baseline_high_gad), bs(1:24,:,baseline_high_gad), [], opt)
legend( {'low baseline gad','model','high baseline gad','model'})
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% model performance in different change in anxiety groups
change_gad = gng_change_over_sessions(rctdata.logtotgad(:,1:4));
change_gad = match_id(change_gad(1,:), ones(length(unique(rctdata.subid)),1), unique(gng.subid), rctdata.subid);
change_gad_split = quantile(change_gad(~ex),2);
got_better = change_gad(~ex) <= change_gad_split(1);
no_change = change_gad(~ex) <= change_gad_split(2);
opt.cl = [0 1 0 0.5]; opt.cl2 = [0 1 0];
figure(); pp.plot_model_performance_basic(as(:,:,got_better), bs(:,:,got_better), [], opt)
opt.cl = [1 0 0 0.5];opt.cl2 = [1 0 0];
hold on; pp.plot_model_performance_basic(as(:,:,no_change), bs(:,:,no_change), [], opt)
legend( {'got better','model','no change','model'})
set(findall(gcf,'-property','FontSize'),'FontSize',18);
