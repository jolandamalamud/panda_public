%% Powersimulations 
number_of_subjects = 172*2; 
d = .1:.1:.7; % true effect size
n_runs = length(d);
pp = power_plot;

%% 1) treatment effect on M or Y
parfor f = 1:1000
    for i = 1:n_runs
        Nsj = number_of_subjects;
        X = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        M = d(i) * X + randn(Nsj,1); % simulated data from LME
        
        lm = fitlm(table(M, X), 'M ~ X');
        est(i,f) = lm.Coefficients.Estimate(2); % fixed effect estimate
        se(i,f) = lm.Coefficients.SE(2); % standard error
        sig(i,f) = double(lm.Coefficients.pValue(2) < 0.05); % significance
    end
end
results.effect = d; results.est = est; results.sig = sig; results.se = se;
figure(); save_and_plot(results, 'drugeffect', 0)
clear results sig se est

%% 2) relation between symptoms Y and mediator M
parfor f = 1:1000
    for i = 1:n_runs
        Nsj = number_of_subjects;
        X = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        M = 0.5 * X + randn(Nsj,1); % mediator
        Y = d(i) * M + 0.5 * X + randn(Nsj,1); % simulated data from LME

        lm = fitlm(table(Y, X, M), 'Y ~ M + X');
        est(i,f) = lm.Coefficients.Estimate(3); % fixed effect estimate
        se(i,f) = lm.Coefficients.SE(3); % standard error
        sig(i,f) = double(lm.Coefficients.pValue(3) < 0.05); % significance
    end
end
results.effect = d; results.est = est; results.sig = sig; results.se = se;
figure(); save_and_plot(results, 'symptomcorrelation', 0)
clear results sig se est

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