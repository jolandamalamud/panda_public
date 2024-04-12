%% load trial data
filepath = '/Users/jolandamalamud/phd/projects/gng_panda_antler/gng_panda/data/';

% Import the data
gngpandadata = readtable([filepath 'gng_panda_data.csv']);

variables = gngpandadata.Properties.VariableNames;

% subject ID
idx_id = find(contains(variables, 'identifier'));
rctdata.subid = table2array(gngpandadata(:,idx_id));

% group allocation
idx_group = find(contains(variables, 'group')); % group allocatio [1=placebo;2=sertraline]
rctdata.group_allocation = table2array(gngpandadata(:,idx_group));

% stratification variables:
idx_cis = find(contains(variables, 'cis')); % CIS-R total score [1=0-11, 2=12-19; 3=>=2-49]
idx_dur = find(contains(variables, 'dur')); % CIS-R depression duration (years) [1=<2years;2=>=2years]
idx_site = find(contains(variables, 'site')); % site [1=Bristol;2=Liverpool;3=York;4=London]
subtable = gngpandadata(:,[idx_cis,idx_dur,idx_site]);
rctdata.strat_variables = table2array(subtable);
rctdata.strat_variables_names = subtable.Properties.VariableNames;

% baseline variables:
idx_marstat = find(contains(variables, 'mar')); % marital status [1=married;2=single;3=separated/divorced/widowed]
idx_age = find(contains(variables, 'age')); % age
idx_deppast = find(contains(variables, 'depressed')); % depression in the past [1=no;2=yes]
idx_edu = find(contains(variables, 'edu')); % highest educational qualification [1=A Level or higher;2=GCSE,standard grade or other;3=no formal qualification]
idx_adpast = find(contains(variables, 'anti')); % antidepressants in the past [1=no;2=yes]
idx_fin = find(contains(variables, 'fin')); % financial difficulty [1=living comfortably or doing alright; 2=just about getting by; 3=finding it difficult or very difficult]
idx_eth = find(contains(variables, 'eth')); % ethinicity [1=white;2:7=ethnic minority]
idx_emp = find(contains(variables, 'emp')); % employment status [1=in paid employment;2=not employed]
idx_sex = find(contains(variables, 'sex')); % sex [1=male;2=female]
subtable = gngpandadata(:, [idx_marstat, idx_age, idx_deppast, idx_edu, ...
    idx_adpast, idx_fin, idx_eth, idx_emp, idx_sex]);
rctdata.baseline_variables = table2array(subtable);
rctdata.baseline_variables_names = subtable.Properties.VariableNames;

% psychiatric scores
week = {'', '_2wk', '_6wk', '_12wk'};
% phq9 log-transformed total scores
idx_totphq = find(strcmp(variables, 'phqtot'));
for t = 2:4
    idx_totphq(t) = find(strcmp(variables, ['phqtot' week{t}]));
end
rctdata.totphq = table2array(gngpandadata(:,idx_totphq)); % PHQ9 total score
rctdata.logtotphq = log10(rctdata.totphq+1); % log-transform
% singe items
rctdata.phqitems = zeros(9,4,size(gngpandadata,1));
for i = 1:9
    for t = 1:4
        idx_phqitem = find(strcmp(variables, ['phq' num2str(i) week{t}]));
        rctdata.phqitems(i, t, :) = table2array(gngpandadata(:,idx_phqitem)); % PHQ9 items
    end
end
rctdata.logphqitems = log10(rctdata.phqitems+1); % log-transform
% anhedonia
rctdata.anhedonia = sq(rctdata.phqitems(1, :, :))'; % PHQ9 item 1 score
rctdata.loganhedonia = sq(rctdata.logphqitems(1, :, :))'; % PHQ9 item 1 score

% gad9 log-transformed total scores
idx_totgad = find(strcmp(variables, 'gadtot'));
for t = 2:4
    idx_totgad(t) = find(strcmp(variables, ['gadtot' week{t}]));
end
rctdata.totgad = table2array(gngpandadata(:,idx_totgad)); % GAD7 total score
rctdata.logtotgad = log10(rctdata.totgad+1);% log-transform
% singe items
rctdata.gaditems = zeros(7,4,size(gngpandadata,1));
for i = 1:7
    for t = 1:4
        idx_gaditem = find(strcmp(variables, ['gad' num2str(i) week{t}]));
        rctdata.gaditems(i, t, :) = table2array(gngpandadata(:,idx_gaditem)); % GAD7 items
    end
end
rctdata.loggaditems = log10(rctdata.gaditems+1); % log-transform

% save 
save([filepath 'gng_preprocessed/rctdata.mat'], 'rctdata');
clear all;