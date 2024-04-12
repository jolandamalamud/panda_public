function data_included = exclude_uninformative_data(data, excluding)
% exclude data due to uninformative task performance
    data_included = data;
    Nsj = size(excluding,1);
    if ndims(data) == 3
        T = size(data,3);
        for t = 1:T
            if ismatrix(data)
                data_included(excluding(:,t) ~= 0, t) = ...
                    nan(sum(excluding(:,t) ~= 0), 1);
            elseif ndims(data) == 3
                data_included(:, excluding(:,t) ~= 0, t) = ...
                    nan(size(data,1), sum(excluding(:,t) ~= 0));
            elseif ndims(data) == 4
                data_included(:, :, excluding(:,t) ~= 0, t) = ...
                    nan(size(data,1), size(data,2), sum(excluding(:,t) ~= 0));
            end
        end
    elseif ismatrix(data)
        if size(data_included,1) == Nsj; data_included = data_included'; end
        data_included(:, excluding ~= 0) = ...
                    nan(size(data,1), sum(excluding ~= 0));
        
    else
        error('unknown data dimensions')
    end
    
end