function errobar_plot(data, c)
    y = sq(nanmean(data,2));
    e = sq(nanstd(data,0,2));
    p = errorbar(y,e,'linewidth', 3, 'color', c); 
    alpha = 0.1;   
    set([p.Bar, p.Line], 'ColorType', 'truecoloralpha', 'ColorData', [p.Line.ColorData(1:3); 255*alpha])
end