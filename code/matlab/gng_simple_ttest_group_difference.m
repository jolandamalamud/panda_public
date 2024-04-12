function [m, t, p] = gng_simple_ttest_group_difference(data, gng, rctdata, opt)

    Nsj = length(rctdata.subid);
    
    if opt.match_data
        data = match_id(data, gng.sess, rctdata.subid, gng.subid);
    end
    if opt.exlude_data
        missing_factor = match_id(gng.ex, gng.sess, rctdata.subid, gng.subid);
        data = exclude_uninformative_data(data, missing_factor ~= 0);
    end
    
    if ndims(data) == 3
        for k = 1:size(data, 3)
            for i = 1:size(data, 1)
                [~,pp,~,stats] = ttest2(data(i,rctdata.group_allocation==2,k), ...
                    data(i,rctdata.group_allocation==1,k));
                p(i,k) = pp; t(i,k) = stats.tstat;
                m(:,i,k) = [nanmean(data(i,rctdata.group_allocation==2,k),2), ...
                	nanmean(data(i,rctdata.group_allocation==1,k),2)];
            end
        end
    elseif ismatrix(data)
        if size(data,1) > 1
            if size(data,1) == Nsj; data = data'; end
            for i = 1:size(data, 1)
                [~,pp,~,stats] = ttest2(data(i,rctdata.group_allocation==2), ...
                    data(i,rctdata.group_allocation==1));
                p(i) = pp; t(i) = stats.tstat;
                m(:,i) = [nanmean(data(i,rctdata.group_allocation==2),2), ...
                	nanmean(data(i,rctdata.group_allocation==1),2)];
            end
        else
            [~,p,~,stats] = ttest2(data(rctdata.group_allocation==2), ...
                data(rctdata.group_allocation==1)); t = stats.tstat;
            m(:,i) = [nanmean(data(rctdata.group_allocation==2)), ...
                nanmean(data(rctdata.group_allocation==1))];
        end
    else
        error('unknown dimenions of data');
    end


end
