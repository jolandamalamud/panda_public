function [datafit] = slope(data)
    
    if length(size(data)) == 3
        Nsj = size(data,2);
        npar = size(data,1);
        datafit = nan(2, npar, Nsj);
    elseif length(size(data)) == 2
        Nsj = size(data,1);
        npar = 1;
        datafit = nan(2,Nsj);
    end
        
    
    for sj = 1:Nsj
        if length(size(data)) == 3
            data_to_fit  = sq(data(:,sj,:));
        elseif length(size(data)) == 2
            data_to_fit = data(sj,:);
        end
        data_to_fit = data_to_fit(:,~isnan(data_to_fit(1,:)));
        if size(data_to_fit,2) > 1
            if length(size(data)) == 3
                for k = 1:npar
                    datafit(:, k, sj) = polyfit([1:size(data_to_fit,2)], ...
                        data_to_fit(k,:), 1)';
            %         y_est(:, k, sj) = polyval(c(:, k, sj),[0:2])';
                end
            elseif length(size(data)) == 2
                datafit(:, sj) = polyfit([1:size(data_to_fit,2)], ...
                        data_to_fit, 1)';
            end
        end
    end

end