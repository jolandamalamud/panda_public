function lme = gng_lme_control_baseline(dependent, independent, rctdata)
    
    [Nsj, T] = size(dependent); 
    
    % transform variables for LME
    subject = repmat(1:Nsj,T,1); subject = subject(:);
    dependent_lme = reshape(dependent', Nsj*T ,1);
    independent_lme = reshape(independent', Nsj*T ,1);
    group = rctdata.group_allocation == 2;
    group_all = [zeros(Nsj,1), repmat(group, 1, T-1)];
    group_lme = reshape(group_all', Nsj*T ,1);
    time = repmat(rctdata.sessions(1:T), 1, Nsj); time = time(:);
    
    control_var = [rctdata.strat_variables, rctdata.baseline_variables];
    control_var_lme = nan(Nsj*T, size(control_var,2));
    for k = 1:size(control_var,2)
        control_var_rep = repmat(control_var(:,k)',T,1); 
        control_var_lme(:, k) = control_var_rep(:);
    end
    
    cistot = control_var_lme(:,1);
    depdur = control_var_lme(:,2);
    site = control_var_lme(:,3);
%     marstat = control_var_lme(:,4); 
    age = control_var_lme(:,5); 
%     deppast = control_var_lme(:,6); 
%     education = control_var_lme(:,7); 
%     AD_past = control_var_lme(:,8);
%     finance = control_var_lme(:,9); 
%     ethnic = control_var_lme(:,10); 
%     empstat = control_var_lme(:,11); 
    sex = control_var_lme(:,12); 

    % run LME
    lme = fitlme(table(dependent_lme, independent_lme, group_lme, ...
        cistot, depdur, site, age, time, sex, subject), ...
        ['dependent_lme ~ independent_lme + group_lme + cistot + ' ...
        'depdur + site + sex + time + age + '...
        '(independent_lme-1|subject) + (1|subject)']);
%     lme = fitlme(table(dependent_lme, independent_lme, group_lme, ...
%         cistot, depdur, site, marstat, age, deppast, education, ...
%         AD_past, finance,  ethnic, empstat, sex, subject), ...
%         ['dependent_lme ~ independent_lme + group_lme + cistot + ' ...
%         'depdur + site + marstat + age + deppast + education + ' ...
%         'AD_past + finance + ethnic + empstat + sex + '...
%         '(independent_lme-1|subject) + (1|subject)']);
end