%% Powersimulations 
number_of_subjects = 172; 
d = .1:.2:.7; % true effect size
n_runs = length(d);
intercept = 0; % fixed intercept
T = 3; % number of sessions

%% 1) treatment effect on Y
% true fixed effect: standardized effect size / pooled SD
% pooled SD = var(intercept) + var(residuals)
a = d * sqrt(2);
parfor f = 1:100
    for i = 1:n_runs
        
        Nsj = number_of_subjects;
        x = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        X = [zeros(Nsj,1), x, x]';
        Y = X .* a(i) + (intercept + randn(1, Nsj)) + randn(3, Nsj); % simulated data from LME
        
        Subject = repmat(1:Nsj,T,1); Subject = Subject(:);
        Y = reshape(Y, Nsj*T ,1);
        X = reshape(X, Nsj*T ,1);
        
        lme = fitlme(table(Y, X, Subject), 'Y ~ X + (1|Subject)');
        a_est(i,f) = lme.Coefficients.Estimate(2); % fixed effect estimate
        sig(i,f) = double(lme.Coefficients.pValue(2) < 0.05); % significance
        ci_a(:,i,f) = [lme.Coefficients.Lower(2),lme.Coefficients.Upper(2)];
        standardized_d(i,f) = lme.Coefficients.Estimate(2) / sqrt(var(lme.randomEffects) + var(lme.residuals));
        
    end
end
results.effect = standardized_d; results.true = a; results.est = a_est; results.sig = sig; results.ci = ci_a;
save_and_plot(results, 'drugeffect', 0)

%% 2) correlation with symptoms Y and mediator M
% true fixed effect: standardized effect size / pooled SD 
% pooled SD = var(intercept) + var(slope) + var(residuals)
b = d * sqrt(3); 
parfor f = 1:1000
    for i = 1:n_runs
       
        Nsj = number_of_subjects;
        x = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        X = [ones(Nsj,1), x, x]';
        M = randn(Nsj,T)'; % mediator: aversive Pavlovian bias
        Y = M .* (b(i) + randn(1, Nsj)) + 0.8 .* X ...
            + (intercept + randn(1, Nsj)) ...
            + randn(3, Nsj); % simulated data from LME
        
        Subject = repmat(1:Nsj,T,1); Subject = Subject(:);
        Y = reshape(Y, Nsj*T ,1);
        X = reshape(X, Nsj*T ,1);
        M = reshape(M, Nsj*T ,1);
        
        % Y^t_n = (i1 + i_n) + (b + b_n) * X^t_n + e1^t_n
        lme = fitlme(table(Y, X, M, Subject), ...
            'Y ~ M + X + (M-1|Subject) + (1|Subject)');
        b_est(i,f) = lme.Coefficients.Estimate(3); % estimated effect between mediator and outcome
        sig(i,f) = double(lme.Coefficients.pValue(3) < 0.05); % significance
        ci_b(:,i,f) = [lme.Coefficients.Lower(3),lme.Coefficients.Upper(3)];
        
    end       
end
results.effect = d; results.true = b; results.est = b_est; results.sig = sig; results.ci = ci_b;
save_and_plot(results, 'symptomcorrelation', 0)

%% 3) mediation power 
alpha = 0.05;
a = d; b = d;
for k = 1:length(d)
    for i = 1:length(d)
        med(k,i) = d(k)*d(i);
        powermed(k,i) = calculate_mediation_power_Sobel(number_of_subjects, a(k), b(i), 1, 1, 1, alpha);
        
    end
end
results.effect = d; results.medeffect = med; results.powermed = powermed;
save_and_plot(results, 'mediation', 1)


function save_and_plot(results, filename, opt)

    save(['matfiles/' filename datestr(now, 'yyyymmdd')], 'results');
    pp = power_plot;
    if opt == 0
        pp.power_plot_lme(results.effect, results.sig);
    elseif opt == 1
        pp.power_plot_med(results.powermed);
    end
    saveas(gcf,['figures/' filename datestr(now, 'yyyymmdd') '.eps']);

end