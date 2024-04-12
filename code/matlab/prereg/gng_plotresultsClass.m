classdef gng_plotresultsClass < matlab.mixin.Copyable
    
    properties
        Ti
        cmap
    end
    
    methods
        
        function obj = gng_plotresultsClass()
            obj.Ti = {'Go to win','Go to avoid','Nogo to win','Nogo to avoid'};
            obj.cmap = colormap(lines);
        end
        
        %% calculate go probability
        function [as, bs] = average_go(obj, data, surr, models)
            as = nan(length(data(1).a)/4,4,length(data));
            if ~isempty(surr)
                bs = nan(length(data(1).a)/4,4,length(data),length(models));
            end
            for sj=1:length(data)
                a = data(sj).a; 
                to_late = data(sj).a ~= 3; 
                a(data(sj).a==3) = nan; a(data(sj).a==2) = 0;
                s = data(sj).s; 
                for ss=1:4
                    i = s==ss; 
                    as(1:sum(i),ss,sj) = a(i);
                end
                if ~isempty(surr)
                    for mdl = 1:length(models)
                        aa = [surr(sj).(models(mdl).name).a];%[surr.(models(mdl).name).surr(sj,:).a]; %[surr(sj).(models(mdl).name).a];
                        b = nan(96,1);
                        b(to_late) = nanmean(aa==1,2);
                        for ss=1:4
                            i = s==ss; 
                            bs(1:sum(i),ss,sj,mdl) = b(i);
                        end
                    end
                end
            end

        end
        
        %% calculate correct probability
        function [pcas, pcbs] = average_correct(obj, data, surr, models)
            cr = [1 1 2 2];
            for sj=1:length(data)
                a = data(sj).a; a(data(sj).a==3) = [];
                s = data(sj).s; s(data(sj).a==3) = [];
                for mdl = 1:length(models)
                    aa = [surr(sj).(models(mdl).name).a];%[surr.(models(mdl).name).surr(sj,:).a]; %[surr(sj).(models(mdl).name).a];
                    for ss=1:4
                        i = s==ss; 
                        pcas(ss,sj) = nanmean(a(i)==cr(ss));
                        pcbs(ss,sj,mdl) = nanmean(aa(i,:)==cr(ss), 'all');
                    end
                end
            end

        end


        function plot_model_performance_basic(obj, as, bs, clust, opt)
            if ~exist('opt') || isempty(opt); 
                cl = 'k'; ll = '-'; lw = 4;
                cl2 = obj.cmap; ll2 = '-'; lw2 = 4;
            else
                cl = opt.cl; ll = opt.ll; lw = opt.lw;
                cl2 = opt.cl2; ll2 = opt.ll2; lw2 = opt.lw2;
            end
            if ~exist('clust') || isempty(clust); c=1; cc=1; clustlabel = '';
            else
                c = clust.c; cc = clust.cc; nc = clust.nc; 
                clustlabel = ['Clustersize = ' num2str(nc)];
            end
            
            for ss=1:4
                subplot(c,4,ss+(cc-1)*4); hon;
                plot(nanmean(as(:,ss, :),3),'color',cl,'linewidth',lw);
                hold on;
                if ~isempty(bs)
                    for mdl = 1:size(bs,4)
                        plot(sq(nanmean(bs(:,ss, :, mdl),3)),'color',cl2(mdl,:),'linestyle', ll2, 'linewidth',lw2);
                    end
                end
                hold on;
                ylim([0 1]);
                title(obj.Ti{ss});
                xlabel('Trial');
                if ss==1; ylabel('Probability Go');end
            end

        end
        
        function plot_model_performance_advanced(obj, as, bs, pcas, pcbs, clust)
            if isempty(clust); c=1; cc=1; clustlabel = '';
            else
                c = clust.c; cc = clust.cc; nc = clust.nc; 
                clustlabel = ['Clustersize = ' num2str(nc)];
            end
                
            nModls = size(bs,4);
            
            subplot(c,5,1+(cc-1)*5)
%             mybar(nanmean(pcas, 2)', obj.cmap(1:4,:));
            mybar(nanmean(pcas, 2)', .7);
            hold on;
            cmp = colormap(lines);
%             er = errorbar(nanmean(pcas,2), nanstd(pcas,0,2));  
%             er.Color = [0 0 0];   
%             er.LineWidth = 3;
%             er.LineStyle = 'none';
            if ~isempty(bs)
            xx = [1:4]'*ones(1,nModls) + ones(4,1)*linspace(-.3,.3,nModls);
            pl = plot(xx, sq(nanmean(pcbs, 2)), '.-', 'markersize', 25, ...
                'linewidth', 3);
            for mdl = 1:nModls
                pl(mdl).Color = [obj.cmap(mdl,:), 0.7];
            end
            end
            hold off;
            xlim([.5 4.5]); ylim([0 1]);
            ylabel({'Probability correct', clustlabel});
            set(gca,'xticklabel',obj.Ti,'xticklabelrotation',30);
            
            for ss=1:4
                subplot(c,5,ss+1+(cc-1)*5);
                p = pcolor(squeeze(as(:,ss,:))'); 
                shading flat;
                colormap([0.75 0.75 0.75; 0.25, 0.25, 0.25])
                hold on;
                plot(nansum(as(:,ss, :),3), 'k','linewidth',5);
                drawnow; set(gca,'YDir','normal');
                yy = [0:size(as,3)/5:size(as,3)]; yy(1) = 1;
                yticks(yy); yticklabels([0:0.2:1]);
                title(obj.Ti{ss});
                hold on;
                if ~isempty(bs)
                for mdl = 1:nModls
                    plot(sq(nansum(bs(:,ss, :, mdl),3)),'linewidth', 5, 'color', [pl(mdl).Color,0.7]);
                end
                end
                title(obj.Ti{ss});
                xlabel('Trial');
                if ss==1; ylabel('Probability Go');end
            end
%             colorbar('Ticks',[0 1], 'TickLabels',{'nogo', 'go'})
    
        end
        
        function plot_empiricalVSgenerated(obj, as, bs)
            
            for ss=1:4
                subplot(2,4,ss);
                imagesc(squeeze(as(:,ss,:))'); 
                newmap1 = contrast(squeeze(as(:,ss,:)));
                colormap(newmap1); hold on;
                plot(nansum(as(:,ss, :),3), 'b','linewidth',4);
                drawnow; set(gca,'YDir','normal');
                yticks(1: (size(as,3)/5)-1:size(as,3)); yticklabels([0:0.2:1]);
                title(obj.Ti{ss});
                if ss==1; ylabel({'Probability Go', 'empirical data'});end
                subplot(2,4,ss+4);
                imagesc(squeeze(bs(:,ss,:))'); 
                plot(nansum(sq(bs(:, ss, :)),2),'linewidth', 4);
                imagesc(squeeze(bs(:,ss,:))'); 
                newmap1 = contrast(squeeze(bs(:,ss,:)));
                colormap(newmap1); hold on;
                plot(nansum(bs(:,ss, :),3), 'b','linewidth',4);
                drawnow; set(gca,'YDir','normal');
                yticks(1: (size(bs,3)/5)-1:size(bs,3)); yticklabels([0:0.2:1]);
                xlabel('Trial');
                if ss==1; ylabel({'Probability Go', 'generated data'});end
            end
            colorbar('Ticks',[0 1], 'TickLabels',{'go', 'nogo'})
    
        end
    end
end