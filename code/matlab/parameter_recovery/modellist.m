model_characteristics(1).name = 'llb';
model_characteristics(1).npar = 1;
model_characteristics(1).labels = {'bias'};

model_characteristics(2).name = 'llba';
model_characteristics(2).npar = 2;
model_characteristics(2).labels = {'\beta','\alpha'};

model_characteristics(3).name = 'llbax';
model_characteristics(3).npar = 3;
model_characteristics(3).labels = {'\beta','\alpha','\gamma'};

model_characteristics(4).name = 'llbaxb';
model_characteristics(4).npar = 4;
model_characteristics(4).labels = {'log \beta','\alpha','\gamma','bias'};

model_characteristics(5).name = 'llbaepxb';
model_characteristics(5).npar = 5;
model_characteristics(5).labels = {'\beta','\alpha','\pi','\gamma','bias'};

model_characteristics(6).name = 'll2baepxb';
model_characteristics(6).npar = 6;
model_characteristics(6).labels = {'\beta_{rew}','\beta_{loss}',...
    '\alpha','\pi','\gamma','bias'};

model_characteristics(7).name = 'll2ba2epxb';
model_characteristics(7).npar = 7;
model_characteristics(7).labels = {'\beta_{rew}','\beta_{loss}',...
'\alpha','\pi_{rew}','\pi_{loss}','\gamma','bias'};

model_characteristics(8).name = 'll2b2a2epxb';
model_characteristics(8).npar = 8;
model_characteristics(8).labels = {'\beta_{rew}','\beta_{loss}',...
    '\alpha_{go}','\alpha_{nogo}','\pi_{rew}','\pi_{loss}','\gamma','bias'};