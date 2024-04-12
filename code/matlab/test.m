filepath = '/Users/jolandamalamud/phd/projects/gng_panda_antler/gng_panda/';
resultsDir = [filepath 'data/modelling_results/sess2'];
D = importdata([filepath 'data/panda_data/panda_datafile_date.mat']);
incl = vertcat(D.sess)=='2';
data = D(incl);

Nsj = length(data);
QQ = nan(2,4,96,Nsj);
VV = nan(4,96,Nsj);
QQs = nan(2,4,24,Nsj);
VVs = nan(4,24,Nsj);
for k = 1:Nsj
    a = data(k).a; 
    r = data(k).r; 
    s = data(k).s; 
    
    VP = zeros(4,1);
    VI = zeros(2,4);
    
    for t = 1:length(a)
        if a(t) == 3; continue; end
        VI(a(t),s(t))   = VI(a(t),s(t)) + (r(t) - VI(a(t),s(t)));
        VP(s(t))        = VP(s(t))      + (r(t) - VP(s(t)));
        QQ(:,:,t,k) = VI; VV(:,t,k) = VP;
    end
    
    for ss = 1:4
        QQs(:,:,:,k) = QQ(:,:,s==ss,k); VVs(:,:,k) = VV(:,s==ss,k);
    end
end

figure()
for ss = 1:4
    subplot(3,4,k); p = plot(sq(nanmean(aac(:,k,:),3))','--','linewidth',3); ylim([0 1]); hold on; 
    subplot(3,4,k+16);p = plot(sq(nanmean(VVs(:,k,:),3))', '--','linewidth',3); 
    if k ==1; ylabel('Q value');end
    subplot(3,4,k+20);p = plot(sq(nanmean(QQs(:,k,:,:),4))', '--','linewidth',3); 
    if k ==1; ylabel('V value');end

end