function lme = run_lme(options, dependent_var, independent_var, labels, ...
    sessions, confounder_tbl)

    Nsj = size(dependent_var, 1);
    T = options.T;
    
    % transform variables for LME
    subject = repmat(1:Nsj,T,1); subject = subject(:);
    time = repmat(sessions(1:T), 1, Nsj); time = time(:);
    
    if size(dependent_var, 2) == T
        dependent_var_lme = reshape(dependent_var', Nsj*T ,1);
    else
        tmp = repmat(dependent_var',T,1); dependent_var_lme = tmp(:);
    end
    
    if exist('confounder_tbl', 'var')
        confounder_tbl_lme = table();
        eq = '';
        for k = 1:size(confounder_tbl,2)
            c = table2array(confounder_tbl(:,k));
            c = repmat(c',T,1); c = c(:);
            confounder_tbl_lme(:,k) = table(c);
            eq = [eq ' + ' confounder_tbl.Properties.VariableNames{k}];
        end
        confounder_tbl_lme.Properties.VariableNames = ...
            confounder_tbl.Properties.VariableNames;
    else
        confounder_tbl_lme = table();
        eq = '';
    end
    
    if length(labels) > 2
        pred_label = labels{2};
        for k = 1:length(labels)-1
            if k < length(labels)-1
                pred_label = [pred_label ' + ' labels{k+2}];
            end
            independent_var_lme(:,k) = reshape(sq(independent_var(:,:,k))', Nsj*T ,1);
        end
    else
        independent_var_lme = reshape(independent_var', Nsj*T ,1);
        pred_label = labels{2};
    end
        
    
    tbl = [array2table([dependent_var_lme, independent_var_lme, time, subject], ...
        'VariableNames',[labels, 'time', 'subject']), confounder_tbl_lme];
    
    if options.random_slopes
        eq = [labels{1} ' ~ ' pred_label eq ' + time + (' labels{2} ...
            '-1|subject) + (1|subject)'];
    elseif options.group_time_interaction
        eq = [labels{1} ' ~ ' pred_label eq ' + time*group + (1|subject)'];
    else
        eq = [labels{1} ' ~ ' pred_label eq ' + time + (1|subject)'];
    end

    % run LME
    lme = fitlme(tbl, eq);

end