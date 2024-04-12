classdef VisualAnalysisPANDA
    %% plots for visual analyses 
    properties
        parameter_label = {'reward_sensitivity', 'loss_sensitivity', ...
            'reward_learning_rate', 'loss_learning_rate', ...
            'appetitiv_Pavlovian_bias', 'aversive_Pavlovian_bias', ...
            'noise', 'go_bias'};
        task_condition_labels = {'go to win', 'go to avoid', ...
            'nogo to win', 'nogo to avoid'};
        alpha = 0.2;
    end
    
    methods
        %% parameter over time
        function data_over_time(obj, data)
            figure();
            for i = 1:8
                subplot(8,1,i); 
                p = errorbar(sq(nanmean(data(i,:,:),2))', ...
                    sq(nanstd(data(i,:,:),0,2))' / ...
                    sqrt(size(data(i,:,:),2)), 'linewidth', 3); 
                set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', ...
                    'ColorData', [p.Line.ColorData(1:3); 255*obj.alpha])
                title(obj.parameter_label(i));
                xlim([0,4]); xticks([]);
            end
            xticks([1:3]); xlabel('sessions')
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
        end
        
        %% psychiatric score over time
        function score_over_time(obj, data, group)
            figure('Position', [0 0 800 600],'defaultTextInterpreter','none');
            for j = 1:2
                gg = j-1; hold on;
                p = errorbar(sq(nanmean(data(group==gg,:),1))', ...
                sq(nanstd(data(group==gg,:),0,1))' / ...
                sqrt(size(data(group==gg,:),1)), 'linewidth', 3); 
                set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', ...
                    'ColorData', [p.Line.ColorData(1:3); 255*obj.alpha])
            end
            xlim([0,5]); xticks([1:4]); xlabel('sessions');
            ylabel('log GAD-7 total score');
            legend( {'placebo','sertraline'})
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
        end
        
        function data_over_time_groupsplit2(obj, data, group, grouplabels, statistics)
            %% parameter over time splitted in groups
            figure('Position', [1922 1 866 1816]); 
            x = [1:size(data,3)]';
            c = ['r', 'b'];
            for i = 1:size(data,1)
                for j = 1:2
                    gg = j-1;
                    subplot(size(data,1),1,i); hold on;
                    y=nan(3,1); e=nan(3,1);
                    if size(group,2) > 1
                        for t = 1:size(data,3)
                            y(t) = sq(nanmean(data(i,group(:,t)==gg,t),2));  
                            e(t) = sq(nanstd(data(i,group(:,t)==gg,t),0,2)/ ... 
                                sqrt(size(data(i,group(:,t)==gg,t),2)));
                        end
                    else
                        y = sq(nanmean(data(i,group==gg,:),2));
                        e = sq(nanstd(data(i,group==gg,:),0,2));
                    end
                    p = errorbar(y,e,'linewidth', 3, 'color', c(j)); 
                    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', ...
                        'ColorData', [p.Line.ColorData(1:3); 255*obj.alpha])
                end
                title(obj.parameter_label(i));
                xlim([0,size(data,3)+1]); xticks([]);
                if exist('statistics', 'var')
                    text(min(xlim)+0.1, mean(ylim), {['r = ' num2str(round(statistics(1,i),2))], ...
                        ['p = ' num2str(round(statistics(2,i),3))]}, 'Horiz','left', 'Vert','bottom');
                end
            end
            xticks([1:3]); xlabel('sessions');
            legend(grouplabels);
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
        end

        function data_over_time_groupsplit4(obj, data, group1, group2, grouplabels)
            figure('Position', [1922 815 866 1002]); 
            x = [1:3]';
            c = {'y','g','r','m'};
            par = [2,4];
            for i = 1:length(par)
                subplot(length(par),1,i);
                for j = 1:2
                    gg = j-1;
                    y=nan(3,1); e=nan(3,1);
                    for t = 1:3
                        y(t) = nanmean(sq(data(par(i),group2(:,t)==0 & ...
                            group1 == gg,t)),2);    
                        e(t) = nanstd(sq(data(par(i),group2(:,t)==0 & ...
                            group1 == gg,t)),0,2) / ...
                            sum(group2(:,t)==0 & group1 == gg);  
                    end
                    p = errorbar(y,e,'linewidth', 3, 'color', c{j});  hold on;  
                    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', ...
                        'ColorData', [p.Line.ColorData(1:3); 255*obj.alpha]);
                    y=nan(3,1); e=nan(3,1);
                    for t = 1:3
                        y(t) = nanmean(sq(data(par(i),group2(:,t)==1 & ...
                            group1 == gg,t)),2);    
                        e(t) = nanstd(sq(data(par(i),group2(:,t)==1 & ...
                            group1 == gg,t)),0,2) / ...
                            sum(group2(:,t)==1 & group1 == gg);  
                    end
                    p = errorbar(y,e,'linewidth', 3, 'color', c{j+2}); 
                    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', ...
                        'ColorData', [p.Line.ColorData(1:3); 255*obj.alpha]);
                end
                    title(obj.parameter_label(par(i)));
                    xlim([0,4]); xticks([]);
            end
            xticks([1:3]); xlabel('sessions');
            legend(grouplabels);
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
        end

        function raw_behaviour_groupsplit2(obj, data, group, missing_factor, grouplabels)
            figure('Position', [1922 815 866 1002]); 
            for t = 1:3
            for k = 1:4
                subplot(3,4,k+(t-1)*4);
                dd = sq(data(:,k,:,t));
                plot(nansum(dd(:,group(:,t)==0)==1,2)/nansum(missing_factor(group(:,t)==0,t)==0),'color', [0.2 0.8 0.2], 'linewidth', 3); hold on; 
                plot(nansum(dd(:,group(:,t)==1)==1,2)/nansum(missing_factor(group(:,t)==1,t)==0), 'color', [0.8 0.2 0.2], 'linewidth', 3); hold on; 
                ylim([0 1]);
                if t == 1
                    title(obj.task_condition_labels{k});
                end
                if k == 1; ylabel({['session ' num2str(t)], 'Go Probability'}); end
            end
            end
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
            xlabel('Trial');
            legend(grouplabels)
        end

        function plot_pca(obj, data, subplot_title)
            % plot eigenvalues loadings
            figure('Position', [-289 901 639 976]); 
            for i = 1:size(data,2)
                subplot(size(data,2),1,i);
                bar(data(:,i)); 
                title(subplot_title(i))
                xticks([]);
            end
            xticks(1:8); xticklabels(obj.parameter_label); xtickangle(45)
            set(findall(gcf,'-property','FontSize'),'FontSize',18);
        end

      
    end
end
