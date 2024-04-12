function [dataTrain, dataTest] = split_train_test(data)
    
    N = size(data,1); 
    % Cross varidation (train: 70%, test: 30%)
    cv = cvpartition(N,'HoldOut',0.3);
    idx = cv.test;
    % Separate to training and test data
    dataTrain = data(~idx,:);
    dataTest  = data(idx,:);

end