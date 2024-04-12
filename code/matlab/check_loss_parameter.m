lower_anx = zeros(Nsj,T);
upper_anx = zeros(Nsj,T);

for t = 1:3
    tertiles = quantile(rctdata.logtotgad(:,t),3);
    lower_anx(:,t) = rctdata.logtotgad(:,t) <= tertiles(1);
    upper_anx(:,t) = rctdata.logtotgad(:,t) >= tertiles(3);
end

loss_p = sq(gng_parameter_included(2,:,:));
loss_lr = sq(gng_parameter_included(4,:,:));

figure(); 
subplot(221); hist(loss_p(lower_anx==1)); xlim([-2,4]);
title('low anxiety'); xlabel('loss sensitivity');
text(min(xlim)+0.1, max(ylim), {['m = ' num2str(round(nanmean(loss_p(lower_anx==1)),2))], ...
        ['std = ' num2str(round(nanstd(loss_p(lower_anx==1)),3))]}, 'Horiz','left', 'Vert','top');
subplot(223); hist(loss_lr(lower_anx==1)); xlabel('loss LR');
text(min(xlim)+0.1, max(ylim), {['m = ' num2str(round(nanmean(loss_lr(lower_anx==1)),2))], ...
        ['std = ' num2str(round(nanstd(loss_lr(lower_anx==1)),3))]}, 'Horiz','left', 'Vert','top');
subplot(222); hist(loss_p(upper_anx==1)); xlabel('loss sensitivity');
title('high anxiety');
text(min(xlim)+0.1, max(ylim), {['m = ' num2str(round(nanmean(loss_p(upper_anx==1)),2))], ...
        ['std = ' num2str(round(nanstd(loss_p(upper_anx==1)),3))]}, 'Horiz','left', 'Vert','top');
subplot(224); hist(loss_lr(upper_anx==1)); xlabel('loss LR');
text(min(xlim)+0.1, max(ylim), {['m = ' num2str(round(nanmean(loss_lr(upper_anx==1)),2))], ...
        ['std = ' num2str(round(nanstd(loss_lr(upper_anx==1)),3))]}, 'Horiz','left', 'Vert','top');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% parameter over time splitted in groups
figure();
x = [1:3]';
c = ['r', 'b'];
for i = 1:3
    if i > 1
        for j = 1:2
            gg = j-1;
            subplot(3,4,1+(i-1)*3); hold on;
            errobar_plot(loss_p(group==gg & lower_anx(:,t)==1,t), ...
                'linewidth', 3, 'color', c(j));
            subplot(3,4,2+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_lr(group==gg & lower_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_lr(group==gg & lower_anx(:,t)==1,t),0,2))),'linewidth', 3, 'color', c(j));
            subplot(3,4,3+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_p(group==gg & upper_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_p(group==gg & upper_anx(:,t)==1,t),0,2))),'linewidth', 3, 'color', c(j));
            subplot(3,4,4+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_lr(group==gg & upper_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_lr(group==gg & upper_anx(:,t)==1,t),0,2))),'linewidth', 3, 'color', c(j));
        end
    else
            subplot(3,4,1+(i-1)*3); hold on;
            errobar_plot(loss_p(lower_anx(:,t)==1,t)', c(1));
            subplot(3,4,2+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_lr(group==gg & lower_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_lr(group==gg & lower_anx(:,t)==1,t),0,2))),'linewidth', 3);
            subplot(3,4,3+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_p(group==gg & upper_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_p(group==gg & upper_anx(:,t)==1,t),0,2))),'linewidth', 3);
            subplot(3,4,4+(i-1)*3); hold on;
            errobar_plot(sq(nanmean(loss_lr(group==gg & upper_anx(:,t)==1,t),2)),...
                sq(nanstd(loss_p(loss_lr(group==gg & upper_anx(:,t)==1,t),0,2))),'linewidth', 3);
    end
%     text(min(xlim)+0.1, mean(ylim), {['r = ' num2str(round(rr_group(1,i),2))], ...
%         ['p = ' num2str(round(rr_group(2,i),3))]}, 'Horiz','left', 'Vert','bottom');
end
xticks([1:3]); xlabel('sessions');
legend('placebo','sertraline');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% Plot behaviour in low and high anxious participants
Ti = {'Go to win','Go to avoid','Nogo to win','Nogo to avoid'};
figure(); 
for t = 1:3
for k = 1:4
    subplot(3,4,k+(t-1)*4);
    dd = sq(raw_behaviour_included(:,k,:,t));
    plot(nansum(dd(:,lower_anx(:,t)==1)==1,2)/nansum(missing_factor(lower_anx(:,t)==1,t)==0),'color', [0.2 0.8 0.2], 'linewidth', 3); hold on; 
    plot(nansum(dd(:,upper_anx(:,t)==1)==1,2)/nansum(missing_factor(upper_anx(:,t)==1,t)==0), 'color', [0.8 0.2 0.2], 'linewidth', 3); hold on; 
    ylim([0 1]);
    if t == 1
        title(Ti{k});
    end
    if k == 1; ylabel({['session ' num2str(t)], 'Go Probability'}); end
end
end
set(findall(gcf,'-property','FontSize'),'FontSize',18);
xlabel('Trial');
legend('low anxiety', 'high anxiety')
