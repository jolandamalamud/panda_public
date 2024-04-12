function [track,l,dl,dsurr] = ll8b2a2epxbdbx(x,D,mu,nui,doprior,options);
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
beta 		= exp(x(1:8));				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(9:10)));		% learning rate
epsilon 	= exp(x(11:12));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
g    		= x(13);                     % irreducible noise 
bias 		= x(14);						% constant bias
deltab      = x(15);
deltax      = x(16);

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
       dQdb = zeros(2,4,8);
       dVdb = zeros(4,8);
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

	rew = (r(t) > 0 & rho == 1) | (r(t) == 0 & rho == 0); % better outcome for specific condition
	er = beta(2*s(t) - (1 * rew)); % outcome sensitivity for each outcome in each condition

	if dodiff
        for k=1:8
			tmp = (dQdb(:,s(t),k) + [epsilon(2-rho)*dVdb(s(t),k);0]);
			dl(k) = dl(k) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
            tmp = 2*s(t) - (1 * rew) == k;
			dQdb(a(t),s(t),k) = (1-alfa(2-rho))*dQdb(a(t),s(t),k) + alfa(2-rho)*er * tmp;
			dVdb(     s(t),k) = (1-alfa(2-rho))*dVdb(     s(t),k) + alfa(2-rho)*er * tmp;
        end

		tmp = (dQde(:,s(t)) + [epsilon(2-rho)*dVde(s(t));0]);
		dl(10-rho) = dl(10-rho) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa(2-rho))*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa(2-rho)*(1-alfa(2-rho));
		dVde(     s(t)) = (1-alfa(2-rho))*dVde(     s(t)) + (er-V(     s(t)))*alfa(2-rho)*(1-alfa(2-rho));

		dl(12-rho) = dl(12-rho) + noise*(p0(a(t))*epsilon(2-rho)*V(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(13) = dl(13) + noise * (1-noise) * (p0(a(t))-1/2)/pg(a(t));
        
        dl(16) = dl(16) + sess(t) * noise * (1-noise) * (p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(14) = dl(14) + noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
        
        dl(15) = dl(15) + sess(t) * noise*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
      
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa(2-rho) * (er - Q(a(t),s(t)));  
	V(s(t))      = V(s(t))      + alfa(2-rho) * (er - V(s(t)     ));
    track.Q(:,:,t) = Q; track.V(:,t) = V;

end

l  = -l; 
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

