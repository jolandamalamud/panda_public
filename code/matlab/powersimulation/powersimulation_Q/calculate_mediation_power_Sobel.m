function power = calculate_mediation_power_Sobel(N, a, b, stdX, stdM, stderr, alpha)

    corrMX = (a * stdX / stdM)^2; % approximated covariance between mediatior and predictor
    stdXerr = stdM * (1 - corrMX); % noise variance in mediator equation
    
    stda = stdXerr / (N * stdX); % standard error of a
    stdb = stderr / (N * stdM * (1-corrMX)); % standard error of b
    
    stdab = sqrt(a^2 * stderr^2 / (N * stdM^2 * (1 - corrMX)) ...
        + b^2 * stdM^2 * (1- corrMX) / (N * stdX^2)); % standard error of a * b
    
    zscore = a * b / stdab;
  
    alpha2 = alpha / 2;
    za2 = norminv(1 - 0.025);
    power = 1 - normcdf(za2 - zscore) + normcdf(- za2 - zscore);
    
end