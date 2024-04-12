function R = gng_parameter_mse_corr(R)
    
    for sj = 1:length(R.Data)
        R.mse(:,sj) = (R.Data(sj).trueParam - R.fit.parest(:,sj)).^2;
        R.scc(:,sj) = corr(R.Data(sj).trueParam, R.fit.parest(:,sj));
        R.trueparams(:,sj) = R.Data(sj).trueParam;
    end

    for p = 1:R.model.npar
        R.pcc(p) = corr(R.trueparams(p,:)', R.fit.parest(p,:)');
    end
    
end