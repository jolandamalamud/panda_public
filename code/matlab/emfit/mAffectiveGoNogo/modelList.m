function model = modelList
% 
% A file like this should be contained in each model class folder and list the
% models to be run, together with some descriptive features. 
% 
% GENERAL INFO: 
% 
% list the models to run here. The models must be defined as likelihood functions in
% the models folder. They must have the form: 
% 
%    [l,dl,dsurr] = ll(parameters,dataToFit,PriorMean,PriorInverseCovariance,doPrior,otherOptions)
% 
% where otherOptions.generatesurrogatedata is a binary flag defining whether to apply the prior, and
% doGenerate is a flag defining whether surrogate data (output in asurr) is
% generated. 
% 
% name: names of model likelihood function in folder models
% npar: number of paramters for each 
% parnames: names of parameters for plotting
% partransform: what to do to the (transformed) parameter estimates to transform them into the
% parameters
%
% Quentin Huys 2018 qhuys@cantab.net

i=0;
% model 1
i=i+1; 
model(i).descr = 'bias model';
model(i).name = 'llb';			
model(i).npar = 1;
model(i).parnames = {'\bias'};
model(i).parnames_untr = {'bias'};
model(i).partransform = {'@(x)'};
% model 2
i=i+1; 
model(i).descr = 'RW model';
model(i).name = 'llba';			
model(i).npar = 2;
model(i).parnames = {'\beta','\alpha'};
model(i).parnames_untr = {'log \beta','siginv \alpha'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))'};
% model 3
i=i+1; 
model(i).descr = 'RW model with irreducible noise';
model(i).name = 'llbax';			
model(i).npar = 3;
model(i).parnames = {'\beta','\alpha','\gamma'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \gamma'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))'};
% model 4
i=i+1; 
model(i).descr = 'RW model with irreducible noise and with constant bias. ';
model(i).name = 'llbaxb';			
model(i).npar = 4;
model(i).parnames = {'\beta','\alpha','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','siginv \gamma','b'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)x'};
% model 5
i=i+1; 
model(i).descr = 'RW model with constant bias towards one action, irreducible noise and positive Pavlovian bias parameter';
model(i).name = 'llbaepxb';			
model(i).npar = 5;
model(i).parnames = {'\beta','\alpha','\pi','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi','siginv \gamma','b'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};
% model 6
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, joint reward/loss sensitivity, and separate positively constrained Pavlovian bias parameters for rewards and losses. ';
model(i).name = 'llba2epxb';			
model(i).npar = 6;
model(i).parnames = {'\beta','\alpha','\pi_{rew}','\pi_{loss}','\gamma','b'};
model(i).parnames_untr = {'log \beta','siginv \alpha','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};
% model 7
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, separate reward and loss sensitivities, and separate positively constrained Pavlovian bias parameters for rewards and losses. ';
model(i).name = 'll2ba2epxb';			
model(i).npar = 7;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha','\pi_{rew}','\pi_{loss}','\gamma','b'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};
% model 8
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxb';			
model(i).npar = 8;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','b'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x'};

nModls = i;

fprintf('%i models in model list\n',i);
