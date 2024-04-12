filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/';

% task data
D = importdata([filepath 'panda_data/panda_datafile_date.mat']);
Nruns = length(D);
T = length(D(1).a);

% rct data
load([filepath 'gng_preprocessed/rctdata.mat']);

subject = [];
choice =  [];
time = [];
stimulus = [];
group = nan(T,Nruns);
cis = nan(T,Nruns);
dep = nan(T,Nruns);
site = nan(T,Nruns);
for sj = 1:Nruns
   subject = [subject; repmat(str2num(D(sj).sjid),T,1)];
   choice =  [choice; D(sj).a];
   time = [time; repmat(str2num(D(sj).sess),T,1)];
   stimulus = [stimulus; D(sj).s];
   idx = find(rctdata.subid == str2num(D(sj).sjid));
   if ~isempty(idx)
       if str2num(D(sj).sess) > 1
           group(:, sj) = rctdata.group_allocation(idx) == 2 * ones(T,1);
       else
           group(:, sj) = zeros(T,1);
       end
       cis(:, sj) = repmat(rctdata.strat_variables(idx,1),T,1);
       dep(:, sj) = repmat(rctdata.strat_variables(idx,2),T,1);
       site(:, sj) = repmat(rctdata.strat_variables(idx,3),T,1);
   end
end

subids = unique(subject);
for sj = 1:length(subids)
    subject(subject == subids(sj)) = sj;
end
choice(choice == 3) = nan; choice(choice == 2) = 0;
required = sum(stimulus == [1,2],2);
valence = sum(stimulus == [1,3],2);
group = group(:);
cis = cis(:);
dep = dep(:);
site = site(:);

tbl = table(subject, choice, time, required, valence, group, cis, dep, site);
    
%% run logistic MLE

glme = fitglme(tbl,['choice ~ 1 + required * valence * time * group + site + cis + dep +' ...
        '(1 + required * valence * time|subject)'], 'Distribution','Binomial');