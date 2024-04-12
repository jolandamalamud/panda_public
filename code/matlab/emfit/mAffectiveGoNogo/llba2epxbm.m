function [l,dl,dsurr] = llba2epxbm(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llba2epxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, constant bias, joint reward/loss sensitivity, and separate positively 
% constrained Pavlovian bias parameters for rewards and losses. 
% 
% Use this within emfit.m to tit RL type models to a group of subjects using EM. 
% 
% Guitart-Masip M, Huys QJM, Fuentemilla L, Dayan P, Duezel E and Dolan RJ
% (2012): Go and nogo learning in reward and punishment: Interactions between
% affect and effect.  Neuroimage 62(1):154-66 
%
% Quentin Huys 2011-2012 qhuys@gatsby.ucl.ac.uk

dodiff=nargout==2;
np = length(x);
beta 		= exp(x(1));				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(2)));		% learning rate
epsilon 	= exp(x(3:4));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
g    		= 1/(1+exp(-x(5)));         % irreducible noise 
bias 		= x(6);						% constant bias 
m           = 1./(1+exp(-x(7)));        % forgetting parameter

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

stimuli = unique(s(~isnan(s)));
actions = unique(a(~isnan(a)));


V=zeros(4,1); 
Q=zeros(2,4); 

dQdb = zeros(2,4);
dQde = zeros(2,4);
dQdm = zeros(2,4);
dVdb = zeros(4,1);
dVde = zeros(4,1);
dVdm = zeros(4,1);


if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
	rho = sum(s(t)==[1 3]);

	q = Q(:,s(t)); 
	q(1) = q(1) + epsilon(2-rho) * V(s(t)) + bias;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
	end
	l = l + log(pg(a(t)));

	er = beta * r(t);

	if dodiff
		tmp = (dQdb(:,s(t)) + [epsilon(2-rho)*dVdb(s(t));0]);
		dl(1) = dl(1) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
		dVdb(     s(t)) = (1-alfa)*dVdb(     s(t)) + alfa*er;
        dQdb(~(actions==a(t)),:) = m*dQdb(~(actions==a(t)),:);
        dQdb(actions==a(t),~(stimuli==s(t))) = m*dQdb(actions==a(t),~(stimuli==s(t)));
        dVdb(~(stimuli==s(t)))= m*dVdb(~(stimuli==s(t)));

		tmp = (dQde(:,s(t)) + [epsilon(2-rho)*dVde(s(t));0]);
		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);
		dVde(     s(t)) = (1-alfa)*dVde(     s(t)) + (er-V(     s(t)))*alfa*(1-alfa);
        dQde(~(actions==a(t)),:) = m*dQde(~(actions==a(t)),:);
        dQde(actions==a(t),~(stimuli==s(t))) = m*dQde(actions==a(t),~(stimuli==s(t)));
        dVde(~(stimuli==s(t)))= m*dVde(~(stimuli==s(t)));

		dl(4-rho) = dl(4-rho) + g*(p0(a(t))*epsilon(2-rho)*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(5) = dl(5) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(6) = dl(6) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        
        tmp = (dQdm(:,s(t)) + [epsilon(2-rho)*dVdm(s(t));0]);
        dl(7) = dl(7) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        dQdm(a(t),s(t)) = (1-alfa)*dQdm(a(t),s(t));   
        dVdm(   s(t)) = (1-alfa)*dVdm(   s(t));   
        dQdm(~(actions==a(t)),:)= m*((1-m)*Q(~(actions==a(t)),:) + dQdm(~(actions==a(t)),:));  
        dQdm(actions==a(t),~(stimuli==s(t))) = m*((1-m)*Q(actions==a(t),~(stimuli==s(t))) + dQdm(actions==a(t),~(stimuli==s(t)))); 
        dVdm(~(stimuli==s(t))) = m*((1-m)*V(~(stimuli==s(t))) + dVdm(~(stimuli==s(t)))); 
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa * (er - V(s(t)     ));
    Q(~(actions==a(t)),:) = m*Q(~(actions==a(t)),:);
    Q(actions==a(t),~(stimuli==s(t))) = m*Q(actions==a(t),~(stimuli==s(t)));
    V(~(stimuli==s(t)))= m*V(~(stimuli==s(t)));
%     track.go(:,t) = q; track.P(:,t) = pg; track.Q(:,:,t) = Q; track.V(:,t) = V;


end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


