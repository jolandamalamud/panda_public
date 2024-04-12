function D_included = exclude_for_fit(D)
    
    cr = [1 1 2 2]; % correct action
   
    for sj = 1:length(D)
        a = D(sj).a;
        as(sj) = nanmean(a(~(a==3)) == 1); 
        for ss = 1:4
            i = D(sj).s == ss;
            pcorr(ss,sj) = nanmean(a(i) == cr(ss));
        end
    end

    % exclude subjects which either only go or nogo
    idx1 = as == 1 | as == 0; % 1=excluded, 0=included
    
    % exclude subjects which perform worse or equal than chance in all conditions
%     idx2 = all(pcorr <= 0.5);

    ex = idx1 | idx2;
    
    D_included = D(~ex);


end
