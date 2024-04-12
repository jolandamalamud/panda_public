function [loadings, v] = dim_reduction_3d(data, dim)

    parcov = cov(reshape(data, size(data,1), size(data,2)*size(data,3))', 'omitrows');
    [v,d] = eig(parcov);
    d = diag(d); d = d(end:-1:1);
    v = v(:,end:-1:1);
    for i = 1:size(data,2)
        for j = 1:size(data,3)
            loadings(:,i,j) = data(:,i,j)' * v(:, 1:dim);
        end
    end
end