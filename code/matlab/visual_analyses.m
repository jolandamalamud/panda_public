%% plots for visual analyses 
parameter_label = {'reward sensitivity', 'loss sensitivity', ...
        'reward learning rate', 'loss learning rate', ...
        'appetitie Pavlovian bias', 'aversive Pavlovian bias', ...
        'noise', 'go bias'};

%% parameter over time

figure(); 
for i = 1:8
    subplot(8,1,i); 
    errorbar(sq(nanmean(parameter(i,:,:),2))', ...
        sq(nanstd(parameter(i,:,:),0,2))', 'linewidth', 3); 
    title(parameter_label(i));
    xlim([0,4]); xticks([]);
end
xticks([1:3]); xlabel('sessions')
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% parameter over time splitted in groups
figure();
x = [1:3]';
c = ['r', 'b'];
for i = 1:8
    for j = 1:2
        gg = j-1;
        subplot(8,1,i); hold on;
        y = sq(nanmean(gng_parameter_included(i,group==gg,:),2));
        e = sq(nanstd(gng_parameter_included(i,group==gg,:),0,2));
        p = errorbar(y,e,'linewidth', 3, 'color', c(j)); 
        alpha = 0.1;   
        set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha])
    end
    title(parameter_label(i));
    xlim([0,4]); xticks([]);
    text(min(xlim)+0.1, mean(ylim), {['r = ' num2str(round(rr_group(1,i),2))], ...
        ['p = ' num2str(round(rr_group(2,i),3))]}, 'Horiz','left', 'Vert','bottom');
end
xticks([1:3]); xlabel('sessions');
legend('placebo','sertraline');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% split in high and low anxiety
figure();
x = [1:3]';
for i = 1:8
    subplot(8,1,i);
    y=nan(3,1); e=nan(3,1);
    for t = 1:3
        y(t) = nanmean(sq(gng_parameter_included(i,lower_anx(:,t)==1,t)),2);    
        e(t) = nanstd(sq(gng_parameter_included(i,lower_anx(:,t)==1,t)),0,2);  
    end
    p = errorbar(y,e,'linewidth', 3, 'color', 'g');  hold on;
    alpha = 0.1;   
    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha])
    y=nan(3,1); e=nan(3,1);
    for t = 1:3
        y(t) = nanmean(sq(gng_parameter_included(i,upper_anx(:,t)==1,t)),2);    
        e(t) = nanstd(sq(gng_parameter_included(i,upper_anx(:,t)==1,t)),0,2);  
    end
    p = errorbar(y,e,'linewidth', 3, 'color', 'r'); 
    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha])
    title(parameter_label(i));
    xlim([0,4]); xticks([]);
end
xticks([1:3]); xlabel('sessions');
legend('lower anxiety','higher anxiety');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% split in high and low anxiety and sertraline and placebo
figure();
x = [1:3]';
c = {'y','g','r','m'};
par = [2,4];
for i = 1:2
    subplot(2,1,i);
    for j = 1:2
        gg = j-1;
        y=nan(3,1); e=nan(3,1);
        for t = 1:3
            y(t) = nanmean(sq(gng_parameter_included(par(i),lower_anx(:,t)==1 & group == gg,t)),2);    
            e(t) = nanstd(sq(gng_parameter_included(par(i),lower_anx(:,t)==1 & group == gg,t)),0,2);  
        end
        p = errorbar(y,e,'linewidth', 3, 'color', c{j});  hold on;
        alpha = 0.2;   
        set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha]);
        y=nan(3,1); e=nan(3,1);
        for t = 1:3
            y(t) = nanmean(sq(gng_parameter_included(par(i),upper_anx(:,t)==1 & group == gg,t)),2);    
            e(t) = nanstd(sq(gng_parameter_included(par(i),upper_anx(:,t)==1 & group == gg,t)),0,2);  
        end
        p = errorbar(y,e,'linewidth', 3, 'color', c{j+2}); 
        set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha]);
    end
        title(parameter_label(par(i)));
        xlim([0,4]); xticks([]);
end
xticks([1:3]); xlabel('sessions');
legend('lower anxiety','lower anxiety 5HT', 'higher anxiety', 'higher anxiety 5HT');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% PCA on cognitive parameters
[pcloadings, pcomponents] = dim_reduction_3d(gng_parameter_included, 1);

figure(); 
for i = 1:3
    subplot(3,1,i);
    bar(pcomponents(:,i)); 
    xticks([]);
end
xticks(1:8); xticklabels(parameter_label); xtickangle(45)
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% split in high and low anxiety and sertraline and placebo
figure();
x = [1:3]';
c = {'g','y','m','r'};
par = [2,4];
for j = 1:2
    gg = j-1;
    y=nan(3,1); e=nan(3,1);
    for t = 1:3
        y(t) = nanmean(sq(pcloadings(1,lower_anx(:,t)==1 & group == gg,t)),2);    
        e(t) = nanstd(sq(pcloadings(1,lower_anx(:,t)==1 & group == gg,t)),0,2);  
    end
    p = errorbar(y,e,'linewidth', 3, 'color', c{j});  hold on;
    alpha = 0.2;   
    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha]);
    y=nan(3,1); e=nan(3,1);
    for t = 1:3
        y(t) = nanmean(sq(pcloadings(1,upper_anx(:,t)==1 & group == gg,t)),2);    
        e(t) = nanstd(sq(pcloadings(1,upper_anx(:,t)==1 & group == gg,t)),0,2);  
    end
    p = errorbar(y,e,'linewidth', 3, 'color', c{j+2}); 
    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha]);
end
    title(parameter_label(par(i)));
    xlim([0,4]); xticks([]);
xticks([1:3]); xlabel('sessions');
legend('lower anxiety','lower anxiety 5HT', 'higher anxiety', 'higher anxiety 5HT');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
