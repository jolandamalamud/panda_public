function T = gng_results_table(tabledata, variablenames, rownames, doprint, tabletitle)
    
    if iscell(tabledata)
        T = table(table(tabledata{1}(:,1),tabledata{2}(:,1),'VariableNames', variablenames{1}), ...
            table(tabledata{1}(:,2),tabledata{2}(:,2),'VariableNames', variablenames{1}), ...
            'VariableNames', variablenames{2}, 'RowNames', rownames);
    elseif ismatrix(tabledata)
        T = array2table(tabledata,'VariableNames', variablenames, ...
                'RowNames', rownames);
    else
        error('wrong data format')
    end
    
    if doprint
        disp('---------------------------------------------------------------')
        if exist('tabletitle'); disp(tabletitle); end
        disp(T)
        disp('---------------------------------------------------------------')
    end

end