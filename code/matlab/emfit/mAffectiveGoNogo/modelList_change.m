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
model(i).descr = 'only bias model, allowing for change in bias over sessions';
model(i).name ='llbdb';			
model(i).npar = 2;
model(i).parnames = {'bias','delta_{bias}'};
model(i).parnames_untr = {'bias', 'delta_{bias}'};
model(i).partransform = {'@(x)x', '@(x)x'};
% model 2
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdepl';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\pi_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\pi_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)exp(x)'};
% model 3
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdepr';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\pi_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\pi_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)exp(x)'};
% model 4
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdb';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x'};
% model 5
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbl';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\beta_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\beta_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)exp(x)'};
% model 6
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbr';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\beta_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\beta_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)exp(x)'};
% model 7
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdal';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\alpha_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\alpha_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)1./(1+exp(-x))'};
% model 8
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdar';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\alpha_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\alpha_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)1./(1+exp(-x))'};
% model 9
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdx';			
model(i).npar = 9;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{\gamma}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{\gamma}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)1./(1+exp(-x))'};
% model 10
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbbr';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\beta_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\beta_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 11
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbbl';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\beta_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\beta_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 12
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbepr';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\pi_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\pi_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 13
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbepl';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\pi_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\pi_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 14
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbar';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\alpha_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\alpha_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 15
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbal';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\pi_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\pi_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 16
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbx';			
model(i).npar = 10;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x'};
% model 17
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxar';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\alpha_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\alpha_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 18
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxal';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\alpha_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\alpha_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 19
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxepr';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\pi_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\pi_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 20
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxepl';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\pi_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\pi_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 21
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxbl';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 22
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxbr';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 23
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxblbr';			
model(i).npar = 12;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\beta_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\beta_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 24
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxblal';			
model(i).npar = 12;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\alpha_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\alpha_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 25
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxblar';			
model(i).npar = 12;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\alpha_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\alpha_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 26
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxblepl';			
model(i).npar = 12;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\pi_{loss}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\pi_{loss}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};
% model 27
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbdbxblepr';			
model(i).npar = 12;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\pi_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\pi_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};
i=i+1; 
model(i).descr = 'RW model with irreducible noise, constant bias, 2 learning rates for reward and loss stimuli, separate reward and loss sensitivities, and separate learning rate and positively constrained Pavlovian bias parameters for rewards and losses.  ';
model(i).name ='ll2b2a2epxbmdbx';			
model(i).npar = 11;
model(i).parnames = {'\beta_{rew}','\beta_{loss}','\alpha_{res}','\alpha_{loss}','\pi_{rew}','\pi_{loss}','\gamma','bias','delta_{bias}','delta_{\gamma}','delta_{\beta_{rew}}','delta_{\pi_{rew}}'};
model(i).parnames_untr = {'log \beta_{rew}','log \beta_{loss}','siginv \alpha_{rew}','siginv \alpha_{loss}','log \pi_{rew}','log \pi_{loss}','siginv \gamma','bias','log delta_{bias}','log delta_{\gamma}','log delta_{\beta_{rew}}','log delta_{\pi_{rew}}'};
model(i).partransform = {'@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)1./(1+exp(-x))','@(x)exp(x)','@(x)exp(x)','@(x)1./(1+exp(-x))','@(x)x','@(x)x','@(x)x','@(x)x','@(x)x'};



nModls = i;

fprintf('%i models in model list\n',i);
