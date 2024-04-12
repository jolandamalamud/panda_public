classdef power_plot
    
   properties
       effect_size
   end
   
   methods

    function power_plot_lme(obj, coeff, sig)

        % plot power as a function of true coefficients
        cmap = colormap(lines); 
        plot(coeff, nanmean(sig,2), 'color', [0,0.3,0.6], 'linewidth', 3); 
        hold on;
        try
            xline(coeff(find(nanmean(sig,2) > 0.95,1)), 'r', 'linewidth', 2)
        end
        yline(0.95, '--', 'color', 'r');
        xlabel('effect size'); 
        ylabel('power'); 
        ylim([0 1]);
        set(findall(gcf,'-property','FontSize'),'FontSize',18);

    end

    function power_plot_med(obj, powermed)

      figure();
      imagesc(powermed); set(gca,'YDir','normal');
      xticks(1:10:length(obj.effect_size)); xticklabels(obj.effect_size(1:10:end)); xlabel('effect between predictor and mediator')
      yticks(1:10:length(obj.effect_size)); yticklabels(obj.effect_size(1:10:end)); ylabel('effect between mediator and outcome')
      c = colorbar;
      c.Label.String = 'Power';
      set(findall(gcf,'-property','FontSize'),'FontSize',18);

    end
    
   end
   
end