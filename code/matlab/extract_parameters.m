%% extract parameters

% load fitted model parameter
change_model = false;
complex_model = 'll2b2a2epxb';
simple_model = 'llb';
filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/';
model.(complex_model) = load([filepath 'modelling_results/final_results/' complex_model '.mat']);
model.(simple_model) = load([filepath 'modelling_results/final_results/' simple_model '.mat']);

% load empirical data
D = importdata([filepath 'panda_data/panda_datafile_date.mat']);
if change_model
    D = combine_sessions(D);
end
bc = gng_basic_characteristics(D);
% D = D(vertcat(D.sess) == '1' | vertcat(D.sess) == '2');
bc = gng_basic_characteristics(D);
% D = D(~(nanmean(bc.pcorr,1) <= 0.5));

% i) exclude subjects where the iL of the most parsimonious model is 3x
% lower than the iL of the random baseline model
idx1 = (model.(complex_model).bf.iL - model.(simple_model).bf.iL) <= 3;
% ii) exclude subjects which either only go or only nogo
for k = 1:length(D)
    if change_model; a = vertcat(D(k).data.a); else; a = D(k).a; end
    aa(k) = nanmean(a(~(a==3)) == 1); 
end
idx2 = aa == 1 | aa == 0;
exlusion_idx = idx1 | idx2; %| nanmean(bc.pcorr,1) < 0.5; % 1=excluded, 0=included
% exlusion_idx = idx1;

if change_model
    gng.subid = str2num(cell2mat({D(:).sjid}'));
    gng.sess = ones(length(gng.subid), 1);
else
    gng.subid = str2num(cell2mat({D(:).sjid}'));
    gng.sess = str2num(cell2mat({D(:).sess}'));
end
gng.params = model.(complex_model).E;
gng.ex = exlusion_idx; % 1=excluded, 0=included
gng.modelfit = model.(complex_model).bf.iL;

in=input(['do you really want to save the extracted parameter from ',complex_model,'? [y for yes, else no]'],'s');
if strcmp(in,'y')
    goalpath = '/Users/jolandamalamud/phd/projects/gng_panda_antler/gng_panda/data/gng_preprocessed/';
    save([goalpath 'gngdata_' complex_model '.mat'], 'gng');
end
