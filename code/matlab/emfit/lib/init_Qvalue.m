function [Q, V, dQdb, dVdb, dQde, dVde, dQdd, dVdd, dQdde, dVdde] ...
    = init_Qvalue()

    V=zeros(4,1); 
    Q=zeros(2,4); 

    dQdb = zeros(2,4);
    dQde = dQdb;
    dQdd = dQdb;
    dQdde = dQdb;
    dVdb = zeros(4,1);
    dVde = dVdb;
    dVdd = dVdb;
    dVdde = dVdb;
    
end