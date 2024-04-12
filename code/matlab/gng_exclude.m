function ex = gng_exclude(simplename, modelname)

    filepath = '~/phd/projects/gng_panda_antler/gng_panda/data/';
    simple = load([filepath 'modelling_results/' simplename '.mat']);
    model = load([filepath 'modelling_results/' modelname '.mat']);
    
    ex = model.bf.iL - simple.bf.iL <= 3;
    
end