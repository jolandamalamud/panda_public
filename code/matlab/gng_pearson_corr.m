function [c, p] = gng_pearson_corr(y, x, gng, rctdata, opt)

    Nsj = length(rctdata.subid);
    
    if opt.match_data
        x = match_id(x, gng.sess, rctdata.subid, gng.subid);
    end
    if opt.exlude_data
        missing_factor = match_id(gng.ex, gng.sess, rctdata.subid, gng.subid);
        x = exclude_uninformative_data(x, missing_factor ~= 0);
    end
    
    if ndims(x) == 3
        for k = 1:size(x, 3)
            for i = 1:size(x, 1)
                [c(:,i,k),p(:,i,k)] = corr(x(i,:,k)', y, 'rows', 'complete');
            end
        end
    elseif ismatrix(x)
        if size(x,1) > 1
            if size(x,1) == Nsj; x = x'; end
            for i = 1:size(x, 1)
                [c(:,i),p(:,i)] = corr(x(i,:)', y, 'rows', 'complete');
            end
        else
            [c,p] = corr(x, y, 'rows', 'complete');
        end
    else
        error('unknown dimenions of data');
    end


end
