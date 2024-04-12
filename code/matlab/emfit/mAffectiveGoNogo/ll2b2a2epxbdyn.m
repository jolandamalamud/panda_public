function [l,dl,dsurr] = ll2b2a2epxbdyn(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = ll2b2a2epxbdyn(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, constant bias, separate reward and loss sensitivities, separate
% learning rate for "go and nogo action", and separate
% positively constrained Pavlovian bias parameters for rewards and losses. 
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
beta 		= exp(x(1:2));				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(3:4)));		% learning rate
epsilon 	= exp(x(5:6));                    % 'pavlovian' parameter. Weigth of Vcue into Qgo
g    		= 1/(1+exp(-x(7)));         % irreducible noise 
bias 		= x(8);						% constant bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 

dQdb = zeros(2,4,2);
dQde = zeros(2,4,2);
dVdb = zeros(4,2);
dVde = zeros(4,2);

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
    rat = sum(s(t) == [1 3]);
	q = Q(:,s(t)); 
	q(1) = q(1) + epsilon(2-rat) * V(s(t)) + bias;    % add Pavlovian effect 
    
	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
	end
	l = l + log(pg(a(t)));
    
    rho = a(t) == 1;
	er = beta(2-rat) * r(t);

	if dodiff
		for k=1:2
			tmp = (dQdb(:,s(t),k) + [epsilon(2-rat)*dVdb(s(t),k);0]);
			dl(k) = dl(k) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
			tmp = (rat==1 & k==1) | (rat==0 & k==2);
			dQdb(a(t),s(t),k) = (1-alfa(2-rho))*dQdb(a(t),s(t),k) + alfa(2-rho)*er * tmp;
			dVdb(     s(t),k) = (1-alfa(2-rho))*dVdb(     s(t),k) + alfa(2-rho)*er * tmp;
            
            tmp = (dQde(:,s(t),k) + [epsilon(2-rat)*dVde(s(t),k);0]);
			dl(k+2) = dl(k+2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
			tmp = (rho==1 & k==1) | (rho==0 & k==2);
            dQde(a(t),s(t),k) = (1-alfa(2-rho))*dQde(a(t),s(t),k) + (er-Q(a(t),s(t)))*alfa(2-rho)*(1-alfa(2-rho))* tmp;
            dVde(     s(t),k) = (1-alfa(2-rho))*dVde(     s(t),k) + (er-V(     s(t)))*alfa(2-rho)*(1-alfa(2-rho))* tmp;
		end
        

		dl(6-rat) = dl(6-rat) + g*(p0(a(t))*epsilon(2-rat)*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(7) = dl(7) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(8) = dl(8) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa(2-rho) * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa(2-rho) * (er - V(s(t)     ));

end
l  = -l; 
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

