function [l,dl,dsurr] = ll2b4a2epxbdbx(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = ll2ba2epxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, constant bias, separate reward and loss sensitivities, and separate
% learning rate and positively constrained Pavlovian bias parameters for rewards and losses. 
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
alfa 		= 1./(1+exp(-x(3:6)));		% learning rate
epsilon 	= exp(x(7:8));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
g    		= x(9);                     % irreducible noise 
bias 		= x(10);						% constant bias
deltab      = x(11);
deltax      = x(12);

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

sess = vertcat(D.data.sess);
a = vertcat(D.data.a); 
r = vertcat(D.data.r); 
s = vertcat(D.data.s); 

rep = length(sess)/length(D.data);
rep = [1, rep, 2*rep];

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
        
    if ismember(t, rep)
       [Q, V, dQdb, dVdb, dQde, dVde] = init_Qvalue();
    end
    
    if isnan(a(t)); continue; end
    
	rho = sum(s(t)==[1 3]);

	q = Q(:,s(t));
    
    q(1) = q(1) + epsilon(2-rho) * V(s(t)) + bias + sess(t) * deltab;
    
	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
    noise = 1/(1+exp(-g-sess(t)*deltax));
	pg = noise*p0 + (1-noise)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t), pemp);
	end
	l = l + log(pg(a(t)));

	er = beta(2-rho) * r(t);

	if dodiff
        
        tmp = (dQdb(:,s(t)) + [epsilon(2-rho)*dVdb(s(t));0]);
        dl(2-rho) = dl(2-rho) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        dQdb(a(t),s(t)) = (1-alfa(s(t)))*dQdb(a(t),s(t)) + alfa(s(t))*er;
        dVdb(     s(t)) = (1-alfa(s(t)))*dVdb(     s(t)) + alfa(s(t))*er;


        tmp = (dQde(:,s(t)) + [epsilon(2-rho)*dVde(s(t));0]);
        dl(s(t)+2) = dl(s(t)+2) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        dQde(a(t),s(t)) = (1-alfa(s(t)))*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa(s(t))*(1-alfa(s(t)));
        dVde(     s(t)) = (1-alfa(s(t)))*dVde(     s(t)) + (er-V(     s(t)))*alfa(s(t))*(1-alfa(s(t)));


		dl(8-rho) = dl(8-rho) + noise*(p0(a(t))*epsilon(2-rho)*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(9) = dl(9) + noise * (1-noise) * (p0(a(t))-1/2)/pg(a(t));
        
        dl(12) = dl(12) + sess(t) * noise * (1-noise) * (p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(10) = dl(10) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        
        dl(11) = dl(11) + sess(t) * noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
      
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa(s(t)) * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa(s(t)) * (er - V(s(t)     ));

end

l  = -l; 
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

