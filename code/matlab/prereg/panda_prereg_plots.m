%% plots for nhb paper - panda gonogo task modelling

% load fitted model parameter
complex_model = 'll2b2a2epxb';
simple_model = 'llb';
filepath = '/Users/jolandamalamud/Jolanda/guitartpanda/';
panda.model.(complex_model) = load([filepath 'prereg_results/panda/' complex_model '.mat']);
panda.model.(simple_model) = load([filepath 'prereg_results/panda/' simple_model '.mat']);

% load empirical data
panda.D = importdata([filepath 'panda_data/panda_datafile_date.mat']);

% load surrogate data
files = dir([filepath 'prereg_results/panda/SurrogateData_*']);

for f= 1:length(files)
    surrogate_data = importdata([filepath 'prereg_results/panda/' files(f).name]);
    model(f).name = files(f).name(15:end-4);

    for s = 1:size(surrogate_data,1)
        surr(s).(model(f).name) = surrogate_data(s,:);
    end
end

% exclude subjects:
idx1 = (panda.model.(complex_model).bf.iL - panda.model.(simple_model).bf.iL <= 3);
% all subjects which either only go or only nogo
for k = 1:length(panda.D)
    a = panda.D(k).a;
    aa(k) = nanmean(a(~(a==3)) == 1); 
end
idx2 = aa == 1 | aa == 0;
idx = idx1 | idx2;

panda.ex_idx = idx;

Ti = {'Go to win', 'Go to avoid', 'Nogo to win', 'Nogo to avoid'};

%% What data to plot?
plot_data = panda.D;
plot_surrdata = surr;
mm = model([3,8,6,1]);

%% BIC plot - changing colors
open([filepath 'prereg_results/panda/figs/bic.fig'])
for s = 1:2
    h = findobj(subplot(1,2,s),'Type','Bar');
    cmp  = colormap(lines);
    h.FaceColor = 'flat';
    h.CData = [cmp(1,:); 0.5 0.5 0.5; 0.5 0.5 0.5; cmp(2,:); cmp(3,:); 0.5 0.5 0.5; 0.5 0.5 0.5; cmp(4,:)];
end
h = findobj(subplot(1,2,1),'Type','Line');
for k = 1:length(h)
    h(k).Color = [h(k).Color, 0.2]; 
end
xlabel('Average model log posterior probability (given prior)');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
f = gcf; f.Position = [200 360 1400 480];
saveas(f, 'prereg_figures/bic', 'png'); savefig(f, 'prereg_figures/bic.fig');

%% Figure to check overall model performance
pp = gng_plotresultsClass; close gcf;
[as, bs] = pp.average_go(plot_data, plot_surrdata, mm);
figure(); pp.plot_model_performance_basic(as, bs, [])
[pcas, pcbs] = pp.average_correct(plot_data, plot_surrdata, mm);
figure(); pp.plot_model_performance_advanced(as, bs, pcas, pcbs,[])

%legend([{''}, {'empirical data'}, {mm.name}]);
sgtitle('');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
f = gcf; f.Position = [200 360 1400 480];
saveas(f, 'prereg_figures/surr', 'png'); savefig(f, 'prereg_figures/surr.fig');

%% Excluding subjects - histogram
figure();
[N,X] = hist(panda.model.ll2b2a2epxb.bf.iL - panda.model.llb.bf.iL,100);
bar(X,N,'facecolor',[0.7 0.7 0.7]); hold on;
xline(3, 'linewidth', 2,'linestyle','--','color', 'r')
ylabel('Participants');
xlabel({'integrated loglikelihood (iL) of most parsimonious model', '- iL of random baseline model'});
set(findall(gcf,'-property','FontSize'),'FontSize',20);
f = gcf; f.Position = [200 360 1400 620];
saveas(f, 'prereg_figures/exclusion_hist', 'png'); savefig(f, 'prereg_figures/exclusion_hist.fig');

%% Excluding subjects - preformance
mm(1).name = complex_model;
mm(2).name = simple_model;
[as, bs] = pp.average_go(plot_data, surr, mm);
figure(); 
for k = 1:4; subplot(1,4,k);
    plot(nanmean(as(:,k,~panda.ex_idx),3),'color', [0.2 0.8 0.2], 'linewidth', 3); hold on; 
    plot(nanmean(as(:,k,panda.ex_idx),3), 'color', [0.8 0.2 0.2], 'linewidth', 3); hold on; 
    plot(sq(nanmean(bs(:,k,~panda.ex_idx,1),3)), '--', 'color', [0.2 0.8 0.2 0.5], 'linewidth', 3); hold on;
    plot(sq(nanmean(bs(:,k,panda.ex_idx,2),3)), '--', 'color', [0.8 0.2 0.2 0.5], 'linewidth', 3); hold on;
    ylim([0 1]);
    title(Ti{k});
    if k == 1; ylabel('Go Probability'); end
    xlabel('Trial');
end
%legend('included','excluded', 'most parsimonious model','random baseline model');
set(findall(gcf,'-property','FontSize'),'FontSize',20);
f = gcf; f.Position = [200 360 1400 620];
saveas(f, 'prereg_figures/exclusion_data', 'png'); savefig(f, 'prereg_figures/exclusion_data.fig');