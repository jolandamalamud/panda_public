%% Powersimulations for registered report NHB - PANDA study

%% Hypothesis 1 - association between sertraline (X) and Pavlovian inhibition (Y)

% Mixed-effect linear regression, regressing aversive Pavlovian bias
% parameter on group allocation (placebo vs sertraline)

number_of_subjects = 439; 
prop_missing = [.3 .32 .36]; % amount of missing data in gng task at different sessions
intercept = 0; % fixed intercept
d = [0:0.01:1]; % true effect size
% true fixed effect: standardized effect size / pooled SD
a = d * sqrt(2);
T = 3; % number of sessions
rng(1);
parfor f = 1:1000
        
    for i = 1:101
        
        Nsj = number_of_subjects;
        x = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        X = [zeros(Nsj,1), x, x]';
        Y = X .* a(i) + (intercept + randn(1, Nsj)) + randn(3, Nsj); % simulated data from LME
        
        % missing data
        for t = 1:T
            Y(t, rand(Nsj,1) < prop_missing(t)) = NaN;
        end
        
        Subject = repmat(1:Nsj,T,1); Subject = Subject(:);
        Y = reshape(Y, Nsj*T ,1);
        X = reshape(X, Nsj*T ,1);
        
        % Y^t_n = (i1 + i_n) + (a + a_n) * X^t_n + e1^t_n (n=subject,t=session)
        lme = fitlme(table(Y, X, Subject), ...
            'Y ~ X + (1|Subject)');
        a_est(i,f) = lme.Coefficients.Estimate(2); % fixed effect estimate
        sig(i,f) = double(lme.Coefficients.pValue(2) < 0.05); % significance
        ci_a(:,i,f) = [lme.Coefficients.Lower(2),lme.Coefficients.Upper(2)];
        
    end
        
end

% save true and estiamted coefficients, CIs and if pvalue < 0.05
save(['matfiles/powersimulations_H1_' datestr(now, 'yyyymmdd')], 'a', 'a_est', 'ci_a', 'sig');

% save figure
pp = power_plot;
pp.power_plot_lme(d, sig);
saveas(gcf,['figures/powersimulations_H1_'  datestr(now, 'yyyymmdd') '.fig']);
saveas(gcf,['figures/powersimulations_H1_'  datestr(now, 'yyyymmdd') '.png']);

clear all, close all;
%% Hypothesis 2 - association between anxiety (Y) and Pavlovian inhibition (M)

% Mixed-effects linear regression regressing psychiatric score on gng task
% parameter controlling for group allocation

number_of_subjects = 439; 
prop_missing = [.3 .32 .36]; % amount of missing data in gng task at different sessions
intercept = 10; % true fixed intercept
d = [0:0.01:1]; % true effect size
% true fixed effect: standardized effect size / pooled SD
b = d * sqrt(3); 
T = 3; % number of sessions
rng(1);

parfor f = 1:1000
    for i = 1:101
       
        Nsj = number_of_subjects;
        x = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        X = [ones(Nsj,1), x, x]';
        M = randn(Nsj,T)'; % mediator: aversive Pavlovian bias
        Y = M .* (b(i) + randn(1, Nsj)) + 0.8 .* X ...
            + (intercept + randn(1, Nsj)) ...
            + randn(3, Nsj); % simulated data from LME
        
        % missing data: 
        for t = 1:T
            M(t, rand(Nsj,1) < prop_missing(t)) = NaN;
        end
        
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
% save true and estimated coefficients, CIs and if pvalue < 0.05
save(['matfiles/powersimulations_H2_' datestr(now, 'yyyymmdd')], 'b', 'b_est', 'ci_b', 'sig');

% save figure
pp = power_plot;
pp.power_plot_lme(d, sig);
saveas(gcf,['figures/powersimulations_H2_' datestr(now, 'yyyymmdd') '.fig']);
saveas(gcf,['figures/powersimulations_H2_' datestr(now, 'yyyymmdd') '.png']);

clear all, close all;

%% Hypothesis 3 - Mediation

Nsj = 439; % number of subjects
alpha = 0.05;
d = [0.01:0.01:1-0.01]; % effect size
a = d; b = d;
for k = 1:length(d)
    for i = 1:length(d)
        med(k,i) = d(k)*d(i);
        powermed(k,i) = calculate_mediation_power_Sobel(Nsj, a(k), b(i), 1, 1, 1, alpha);
        
    end
end
% save power for mediation effect
save(['matfiles/powersimulations_H2_' datestr(now, 'yyyymmdd')], 'powermed', 'a', 'b', 'med');

% save figure
pp = power_plot;
pp.effect_size = d;
pp.power_plot_med(powermed);
saveas(gcf,['figures/powersimulations_H3_' datestr(now, 'yyyymmdd') '.fig']);
saveas(gcf,['figures/powersimulations_H3_' datestr(now, 'yyyymmdd') '.png']);

clear all, close all;
%% Hypothesis 4 - association between anxiety (Y) and Pavlovian inhibition (M)

number_of_subjects = 439; 
prop_missing = [.3 .32 .36]; % amount of missing data in gng task at different sessions
intercept = 10; % true fixed intercept
c = [0:0.01:1]; % true effect size
T = 3; % number of sessions
rng(1);

parfor f = 1:1000
    for i = 1:101
       
        Nsj = number_of_subjects;
        X = double(rand(Nsj, 1) > 0.5); % predictor: group allocation
        M = randn(Nsj, 1); % mediator: aversive Pavlovian bias
        Y = 0.16 * M + 0.8 * X + c(i) * X.*M + randn(Nsj,1); % simulated data from LME
        
        % missing data: 
        M(rand(Nsj,1) < prop_missing(1)) = NaN;
        
        lm = fitlm(table(Y, X, M), 'Y ~ M + X + X * M');
        c_est(i,f) = lm.Coefficients.Estimate(4); % estimated effect between mediator and outcome
        sig(i,f) = double(lm.Coefficients.pValue(4) < 0.05); % significance
        
    end
        
end
% save true and estimated coefficients, CIs and if pvalue < 0.05
save(['matfiles/powersimulations_H4_' datestr(now, 'yyyymmdd')], 'c', 'c_est', 'sig');

% save figure
pp = power_plot;
pp.power_plot_lme(c, sig);
saveas(gcf,['figures/powersimulations_H4_' datestr(now, 'yyyymmdd') '.fig']);
saveas(gcf,['figures/powersimulations_H4_' datestr(now, 'yyyymmdd') '.png']);

clear all, close all;