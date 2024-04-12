function gng_plot_parameter_recovery(data, model, plot1)
    figure();
    % plot parameter recovery for different models
    if plot1
        for k = 1:size(data.partrue,1)
            npar(k) = length(data.partrue(data.partrue(:,1,k) ~= 0,1,k));
        end
        [~,I] = sort(npar, 'ascend');
        for i = 1:length(model)
            ii = find([model.npar] == npar(I(i)));
            partrue = data.partrue(data.partrue(:,1,I(i)) ~= 0,:,I(i));
            parest = data.parest(data.parest(:,1,I(i)) ~= 0,:,I(i));
            for p = 1:npar(I(i))
                subplot(length(model),max(npar),p+(i-1)*8); 
                scatter(partrue(p,:), parest(p,:));
                xlabel(['true ' model(ii).labels{p}]);
                ylabel(['estimated ' model(ii).labels{p}]);
            end
        end
        disp('parameter correlations:'); disp(data.pcc(I,:));
    % plot parameter recovery for winning models
    else
        bar([1:size(data,1)],nanmean(data,2));              
        hold on
        er = errorbar(nanmean(data,2), nanstd(data,0,2));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none'; 
        title('winning model');
        ylabel({'correlation between true', 'and estimated parameters'});
        set(gca, 'FontSize', 18, 'xticklabel', model.labels, 'xticklabelRotation', 45);
    end
end