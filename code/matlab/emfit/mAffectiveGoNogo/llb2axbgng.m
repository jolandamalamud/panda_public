function [l,dl,dsurr] = llb2axbgng(x,D,mu,nui,doprior,options);
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
alfa 		= x(2);		% learning rate
g    		= 1/(1+exp(-x(3)));         % irreducible noise 
bias 		= x(4);						% constant bias 
kappa       = x(5);

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

V=zeros(1,4); 
Q=zeros(2,4); 

dQdb = zeros(2,4);
dQde = zeros(2,4);
dQda = zeros(2,4);
dVdb = zeros(4,1);
dVde = zeros(4,1);
dVda = zeros(4,1);

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	dodiff=0;
end

for t=1:length(a)
	rho = sum(s(t)==[1 3]);

	q = Q(:,s(t)); 
	q(1) = q(1) + bias;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
	end
	l = l + log(pg(a(t)));
    
    if r(t) == 1 && a(t) == 1; rat = 1;
    elseif r(t) == -1 && a(t) == 2; rat = -1;
    else rat = 0; 
    end
        
    lr = 1./(1+exp(- alfa - rat * kappa));
	er = beta * r(t);

	if dodiff
		tmp = dQdb(:,s(t));
		dl(1) = dl(1) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQdb(a(t),s(t)) = (1-lr)*dQdb(a(t),s(t)) + lr*er;
		dVdb(     s(t)) = (1-lr)*dVdb(     s(t)) + lr*er;

		tmp = dQde(:,s(t));
		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-lr)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*lr*(1-lr);
		dVde(     s(t)) = (1-lr)*dVde(     s(t)) + (er-V(     s(t)))*lr*(1-lr);

		dl(3) = dl(3) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(4) = dl(4) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        
        tmp = dQda(:,s(t));
		dl(5) = dl(5) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQda(a(t),s(t)) = (1-lr)*dQda(a(t),s(t)) + (er-Q(a(t),s(t)))*rat*lr*(1-lr);
		dVda(     s(t)) = (1-lr)*dVda(     s(t)) + (er-V(     s(t)))*rat*lr*(1-lr);
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + lr * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + lr * (er - V(s(t)     ));
    
%     track.go(:,t) = q; track.P(:,t) = pg; track.Q(:,:,t) = Q; track.V(:,t) = V;


end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


