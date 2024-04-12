function tab = ttest_missingness(data, missing_factor)

    T = size(missing_factor,2);
    Np = size(data,1);
    
    tab = table();
    for t = 1:size(missing_factor,2)
        if length(size(data)) == 3
            d1 = sq(data(missing_factor(:,t)==0,t,:));
            d2 = sq(data(missing_factor(:,t)==1,t,:));
        elseif length(size(data)) == 2
            d1 = data(missing_factor(:,t)==0,:);
            d2 = data(missing_factor(:,t)==1,:);
        end
        m1 = sq(nanmean(d1,1));
        m2 = sq(nanmean(d2,1));
        [h, pval] = ttest2(d1, d2, 0.05/(Np*3));
        tab = [tab, table(h', pval', m1', m2', 'VariableNames', ...
            {['sig T' num2str(t)], ['pval T' num2str(t)], ['not miss T' num2str(t)], ['miss T' num2str(t)]})];
    end
    
end