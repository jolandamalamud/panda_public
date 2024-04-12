function estimate = retrieve_estimate_lme(tbl, lme_function)

    lme = fitlme(tbl, lme_function);
    estimate = lme.Coefficients.Estimate;

end