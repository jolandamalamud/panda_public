% investigate mean Q and V values
filepath = '/Users/jolandamalamud/phd/projects/gng_panda_antler/gng_panda/';
resultsDir = [filepath 'data/modelling_results/sess2'];
D = importdata([filepath 'data/panda_data/panda_datafile_date.mat']);
incl = vertcat(D.sess)=='2';
data = D(incl);

% % deal with trials where participants pressed to late
if any([data.a] == 3, 'all')
    for k = 1:length(data)
%         data(k).a(data(k).a == 3) = 1; % assume them as go trials
        to_late = data(k).a == 3; % exclude those trials
        data(k).a = data(k).a(~to_late);
        data(k).s = data(k).s(~to_late);
        data(k).r = data(k).r(~to_late);
        data(k).Nch = length(data(k).a);
    end
end

model = 'll2b2a2epxb';
load([resultsDir '/' model '.mat'])
fstr = str2func(model);
T = 96;

doprior = 0;
options.generatesurrogatedata=0;
for k = 1:length(data)
    QV_data(k) = fstr([E(:,k)],data(k),[],[],doprior,options);
end
%%
nsamples = 10;
options.generatesurrogatedata=1;
for i = 1:nsamples
    for k = 1:length(data)
        [QV_surr(k,i),~,~,surr(k,i)] = fstr([E(:,k)],data(k),[],[],doprior,options);
    end
end

%%
data = D(incl);
QQ = nan(2,4,T,length(data),2);
CC = nan(2,4,T,length(data),2); 
VV = nan(4,T,length(data),2);
PP = nan(2,T,length(data),2);
GG = nan(2,T,length(data),2);
% LL = nan(4,T,length(data),2);
% WW = nan(T,length(data),2);
for k = 1:length(data)
    i = data(k).a ~= 3; 
    QQ(:,:,i,k,1) = QV_data(k).Q; 
%     CC(:,:,i,k,1) = QV_data(k).C; 
    VV(:,i,k,1) = QV_data(k).V; 
    PP(:,i,k,1) = QV_data(k).P;
    GG(:,i,k,1) = QV_data(k).go;
%     LL(:,i,k,1) = QV_data(k).L;
%     WW(i,k,1) = QV_data(k).w;
    for s = 1:nsamples
        qq(:,:,:,s) = QV_surr(k,s).Q;
%         cc(:,:,:,s) = QV_surr(k,s).C;
        vv(:,:,s) = QV_surr(k,s).V;
        pp(:,:,s) = QV_surr(k,s).P;
        gg(:,:,s) = QV_surr(k,s).go;
%         ll(:,:,s) = QV_surr(k,s).L;
%         ww(:,s) = QV_surr(k,s).w;
    end
        
    QQ(:,:,i,k,2) = nanmean(qq,4);
    VV(:,i,k,2) = nanmean(vv,3);
%     CC(:,:,i,k,2) = nanmean(cc,4);
    PP(:,i,k,2) = nanmean(pp,3);
    GG(:,i,k,2) = nanmean(gg,3);
%     LL(:,i,k,2) = nanmean(ll,3);
%     WW(i,k,2) = nanmean(ww,2);
    clear qq vv pp gg ll ww cc
end
%%
aa = nan(T,length(data)); 
aas = aa;
for sj = 1:length(data)
    i = data(sj).a ~= 3; 
    aa(i,sj) = data(sj).a(i)==1; 
    aas(i,sj) = nanmean([surr(sj,:).a]==1,2); 
end

aac = nan(24,4,length(data)); 
aasc = aac;
gg = nan(2,24,4,length(data),2); 
pg = nan(2,24,4,length(data),2); 
% ll = nan(4,24,4,length(data),2); 
% ww = nan(24,4,length(data),2); 
for sj = 1:length(data) 
    s = data(sj).s; 
    for ss = 1:4 
        i = s == ss; 
        aac(:,ss,sj) = aa(i,sj); 
        aasc(:,ss,sj) = aas(i,sj);
        gg(:,:,ss,sj,:) = GG(:,i,sj,:);
        pg(:,:,ss,sj,:) = PP(:,i,sj,:);
%         ll(:,:,ss,sj,:) = LL(:,i,sj,:);
%         ww(:,ss,sj,:) = WW(i,sj,:);
    end
end
%%
% vap = VisualAnalysisPANDA;
% 
% figure(); 
% for k = 1:4
%     subplot(5,4,k); p = plot(sq(nanmean(aac(:,k,:),3))','--','linewidth',3); hold on; 
%     plot(sq(nanmean(aasc(:,k,:),3))','linewidth',3,'color',[p.Color, 0.5]); ylim([0,1]);
%     if k ==1; ylabel('go choice'); legend('empirical', 'surrogate');end
%     subplot(5,4,k+4);p = plot(sq(nanmean(ww(k,:,1),3))', '--','linewidth',3); 
%     hold on; plot(sq(nanmean(ww(:,k,:,2),3))','linewidth',3,'color',[p.Color, 0.5]);
%     if k ==1; ylabel('Pavlovian weight'); legend('go', 'nogo');end
%     subplot(5,4,k+8);p = plot(sq(nanmean(ll(k,:,k,:,1),4))', '--','linewidth',3); 
%     hold on; plot(sq(nanmean(ll(k,:,k,:,2),4))','linewidth',3,'color',[p.Color, 0.5]);
%     if k ==1; ylabel('go weight');end
%     subplot(5,4,k+12);p = plot(sq(nanmean(QQ(:,k,:,:,1),4))', '--','linewidth',3); 
%     hold on; for i = 1:2; plot(sq(nanmean(QQ(i,k,:,:,2),4))','linewidth',3,'color',[p(i).Color, 0.5]); end
%     if k ==1; ylabel('Q value');end
%     subplot(5,4,k+16);p = plot(sq(nanmean(VV(k,:,:,1),3))', '--','linewidth',3); 
%     hold on; plot(sq(nanmean(VV(k,:,:,2),3))','linewidth',3,'color',[p.Color, 0.5]);
%     if k ==1; ylabel('V value');end
% end
% sgtitle(model);
% set(findall(gcf,'-property','FontSize'),'FontSize',18);
%%
vap = VisualAnalysisPANDA;

figure(); 
for k = 1:4
    subplot(5,4,k); p = plot(sq(nanmean(aac(:,k,:),3))','--','linewidth',3); hold on; 
    plot(sq(nanmean(aasc(:,k,:),3))','linewidth',3,'color',[p.Color, 0.5]); ylim([0,1]);
    if k ==1; ylabel('go choice'); legend('empirical', 'surrogate');end
    subplot(5,4,k+4);p = plot(sq(nanmean(pg(:,:,k,:,1),4))', '--','linewidth',3); 
    hold on; for i = 1:2; plot(sq(nanmean(pg(i,:,k,:,2),4))','linewidth',3,'color',[p(i).Color, 0.5]); end
    if k ==1; ylabel('go probability'); legend('go', 'nogo');end
    subplot(5,4,k+8);p = plot(sq(nanmean(gg(:,:,k,:,1),4))', '--','linewidth',3); 
    hold on; for i = 1:2; plot(sq(nanmean(gg(i,:,k,:,2),4))','linewidth',3,'color',[p(i).Color, 0.5]); end
    if k ==1; ylabel('go weight');end
    subplot(5,4,k+12);p = plot(sq(nanmean(QQ(:,k,:,:,1),4))', '--','linewidth',3); 
    hold on; for i = 1:2; plot(sq(nanmean(QQ(i,k,:,:,2),4))','linewidth',3,'color',[p(i).Color, 0.5]); end
    if k ==1; ylabel('Q value');end
    subplot(5,4,k+16);p = plot(sq(nanmean(VV(k,:,:,1),3))', '--','linewidth',3); 
    hold on; plot(sq(nanmean(VV(k,:,:,2),3))','linewidth',3,'color',[p.Color, 0.5]);
    if k ==1; ylabel('V value');end
end
sgtitle(model);
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% 
figure(); 
for k =1:4; subplot(1,4,k);p = plot(sq(nanmean(aac(:,k,:),3))','--','linewidth',3); hold on;
plot(sq(nanmean(pg(1,:,k,:,1),4))','linewidth',3,'color',[p.Color, 0.5]); ylim([0,1]); end
%%
subs = [1:3];
figure(); 
for sj = 1:length(subs)
for k = 1:4; subplot(length(subs),4,k+(sj-1)*4);
    bar(sq(aac(:,k,subs(sj)))'); hold on;
%     plot(sq(aasc(:,k,subs(sj)))','linewidth',3); ylim([0,1]);
    plot(sq(pg(1,:,k,subs(sj),1))','linewidth',3); ylim([0,1]); 
end
end
