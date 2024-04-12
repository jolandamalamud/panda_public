function results = lme_fu_separate(dependent, independent, options)
% LME for baseline and one of the two follow-up separate
    Nsj = size(dependent,1);
    group_over_time = [zeros(Nsj,1), repmat(options.group, 1, options.T-1)];

    for k = 1:2
        disp(['LME including baseline and follow-up ' num2str(k) ':'])
        if isempty(independent)
            independent_variable = group_over_time;
        else
            independent_variable = ...
                cat(3,independent(:,[1,k+1]), group_over_time);
        end
        % run LME
        lme = run_lme(options, dependent(:,[1,k+1]), ...
            independent_variable, options.labels, ...
            options.time([1, k+1]), options.confounders);
        results{k} = lme;
        if options.disp; disp(lme); end
        disp('--------------------------------------------------------')
    end

end