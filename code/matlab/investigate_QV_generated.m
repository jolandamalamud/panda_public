model = 'll2aepxbcon';
fstr = str2func(model);
T = 96;
Nsj = 100;

doprior = 0;
options.generatesurrogatedata=1;

setparams.mu = [2,2,0.2,0,1]; setparams.sigma = 0.001 * ones(length(setparams.mu));
E =  mvnrnd(setparams.mu, setparams.sigma, Nsj)';

for k = 1:Nsj
    data(k).a = zeros(1,T);
	data(k).r = zeros(1,T);
	rs = randperm(T);	
	s = [1:4]'*ones(1,T/4);
    data(k).s = s(rs);		
    [QV_data(k),~,~,surr(k)] = fstr([E(:,k)],data(k),[],[],doprior,options);
end
%%
VI = nan(2,4,T,Nsj,2); 
VQ = nan(4,T,Nsj,2);
LL = nan(4,T,Nsj,2);
WW = nan(T,Nsj,2);
UU = nan(T,Nsj,2);
for k = 1:Nsj
    i = data(k).a ~= 3; 
    VI(:,:,i,k,1) = QV_data(k).Q; 
    VP(:,i,k,1) = QV_data(k).V; 
    LL(:,i,k,1) = QV_data(k).L;
    WW(i,k,1) = QV_data(k).w;
    UU(i,k,1) = QV_data(k).u;
end
%%
aa = nan(T,Nsj); 
aas = aa;
for sj = 1:Nsj
    i = data(sj).a ~= 3; 
    aa(i,sj) = data(sj).a(i)==1; 
end

aac = nan(24,4,Nsj); 
ww = nan(24,4,Nsj,2); 
uu = nan(24,4,Nsj,2); 
for sj = 1:Nsj
    s = data(sj).s; 
    for ss = 1:4 
        i = s == ss; 
        aac(:,ss,sj) = aa(i,sj); 
        ww(:,ss,sj,:) = WW(i,sj,:);
        uu(:,ss,sj,:) = UU(i,sj,:);
    end
end
%%
vap = VisualAnalysisPANDA;

figure(); 
for k = 1:4
    subplot(6,4,k); p = plot(sq(nanmean(aac(:,k,:),3))','--','linewidth',3); ylim([0 1]); hold on; 
    subplot(6,4,k+4);p = plot(sq(nanmean(uu(:,k,:,1),3))', '--','linewidth',3); 
    if k ==1; ylabel('update');end
    subplot(6,4,k+8);p = plot(sq(nanmean(LL(k,:,:,1),3))', '--','linewidth',3); 
    if k ==1; ylabel('L');end
    subplot(6,4,k+12);p = plot(sq(nanmean(ww(:,k,:,1),3))', '--','linewidth',3); 
    if k ==1; ylabel('w');end
    subplot(6,4,k+16);p = plot(sq(nanmean(VI(:,k,:,:,1),4))', '--','linewidth',3); 
    if k ==1; ylabel('Q value');end
    subplot(6,4,k+20);p = plot(sq(nanmean(VP(k,:,:,1),3))', '--','linewidth',3); 
    if k ==1; ylabel('V value');end
end
sgtitle(model);
set(findall(gcf,'-property','FontSize'),'FontSize',18);