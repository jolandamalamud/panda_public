%% mediation analysis using boostrapping

% Y = c*X + er

% M = a*X + er
options.T = 2;
options.random_slopes = false;
options.interaction = false;
loss_lR = sq(gng_parameter_included(4,:,1:options.T));
group_over_time = [zeros(Nsj,1), repmat(group, 1, options.T-1)];

bootstat1 = bootstrap_lme(options, loss_lR, group_over_time, ...
    {'LRloss', 'group'}, [0,2,6], stratification_tbl);

stats.H3.a = mean(bootstat1.estimates(:,2));
stats.H3.ci_a = calculate_CI(bootstat1.estimates(:,2), 95, true);

% Y = b*M + c*X + th*M*X + er
options.T = 3;
options.interaction = true;
interaction = 'LRloss*group';
loss_lR = sq(gng_parameter_included(4,:,1:options.T));
group_over_time = [zeros(Nsj,1), repmat(group, 1, options.T-1)];
independent_variables = cat(3, loss_lR, group_over_time);

bootstat2 = bootstrap_lme(options, rctdata.logtotgad(:,1:options.T), ...
    independent_variables, {'anxiety', 'LRloss', 'group'}, [0,2,6], ...
    [stratification_tbl, baseline_variables_tbl], interaction);

stats.H3.b = mean(bootstat2.estimates(:,2));
stats.H3.ci_b = calculate_CI(bootstat2.estimates(:,2), 95, true);
stats.H3.th = mean(bootstat2.estimates(:,end));
stats.H3.ci_th = calculate_CI(bootstat2.estimates(:,end), 95, true);

% pure natural indirect effect -> in control group
stats.H3.pnie = mean(bootstat1.estimates(:,2) .* bootstat2.estimates(:,2));
stats.H3.ci_ab_pnie = calculate_CI(bootstat1.estimates(:,2) .* ...
bootstat2.estimates(:,2), 95, true); % 95% confidence interval

% total natural indirect effect -> in treatment group
stats.H3.tnie = mean(bootstat1.estimates(:,2) .* bootstat2.estimates(:,2) + ...
    bootstat1.estimates(:,2) .* bootstat2.estimates(:,end));
stats.H3.ci_ab_tnie = calculate_CI(bootstat1.estimates(:,2) .* bootstat2.estimates(:,2) + ...
    bootstat1.estimates(:,2) .* bootstat2.estimates(:,end), 95, true); % 95% confidence interval
