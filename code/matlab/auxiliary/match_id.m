function [data_matched] = match_id(data, sessions, id1, id2)
% match RCT ids (id1) with cognitive task data ids (id2)
    T = max(sessions); 
    Nsj = length(id1);
    if ndims(data) == 3
        data_matched = nan(size(data,1), size(data,2), Nsj, T);
    elseif ndims(data) == 2
        if all(size(data) > 1)
            data_matched = nan(size(data,1), Nsj, T);
        elseif ndims(data) == 2
            data_matched = nan(Nsj, T);
        end
    else 
        error('wrong data type');
    end
    
    for s = 1:Nsj
        subsessions = sessions(id2 == id1(s));
        if ~isempty(subsessions)
            if ndims(data) == 3
                data_matched(:, :, s, subsessions) = data(:,:,id2 == id1(s));
            elseif ndims(data) == 2
                if all(size(data) > 1)
                    data_matched(:, s, subsessions) = data(:,id2 == id1(s));
                elseif ndims(data) == 2
                    data_matched(s, subsessions) = data(id2 == id1(s));
                end
            end
        end
    end
    
end