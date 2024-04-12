function lme = gng_lme_control_strat(dependent, rctdata)
    
    [Nsj, T] = size(dependent); 
    
    % transform variables for LME
    subject = repmat(1:Nsj,T,1); subject = subject(:);
    dependent_lme = reshape(dependent', Nsj*T ,1);
    group = rctdata.group_allocation == 2;
    group_all = [zeros(Nsj,1), repmat(group, 1, T-1)];
    group_lme = reshape(group_all', Nsj*T ,1);
    time = repmat(rctdata.sessions(1:T), 1, Nsj); time = time(:);
    
    control_var = rctdata.strat_variables;
    control_var_lme = nan(Nsj*T, size(control_var,2));
    for k = 1:size(control_var,2)
        control_var_rep = repmat(control_var(:,k)',T,1); 
        control_var_lme(:, k) = control_var_rep(:);
    end
    
    cistot = control_var_lme(:,1);
    depdur = control_var_lme(:,2);
    site = control_var_lme(:,3);

    % run LME
    lme = fitlme(table(dependent_lme, group_lme, time, ...
        cistot, depdur, site, subject), ...
        'dependent_lme ~ group_lme + time + cistot + depdur + site + (1|subject)');
end