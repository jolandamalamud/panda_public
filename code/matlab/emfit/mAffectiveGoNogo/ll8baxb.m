function [l,dl,dsurr] = ll8baxb(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = ll2baxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with irreducible
% noise, constant bias and separate reward and loss sensitivities
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
beta 		= x(1:8);				% sensitivity to rewards and losses
alfa 		= 1./(1+exp(-x(9)));		% learning rate
g    		= 1/(1+exp(-x(10)));		% irreducible noise 
bias 		= x(11);						% constant bias 

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
        V=zeros(1,4); 
        Q=zeros(2,4); 

        dQdb = zeros(2,4,8);
        dQde = zeros(2,8);
        dVdb = zeros(4,4);
        dVde = zeros(4,1);
    end
    
    if isnan(a(t)); continue; end
	rho = sum(s(t)==[1 3]); % rewarding stimuli

	q = Q(:,s(t)); 
	q(1) = q(1)  + bias;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t), pemp);
	end
	l = l + log(pg(a(t)));
    
	rew = (r(t) > 0 & rho == 1) | (r(t) == 0 & rho == 0); % better outcome for specific condition
	er = beta(2*s(t) - (1 * rew)); % outcome sensitivity for each outcome in each condition

	if dodiff
        for k=1:8
			tmp = dQdb(:,s(t),k);
			dl(k) = dl(k) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
			tmp = 2*s(t) - (1 * rew) == k;
			dQdb(a(t),s(t),k) = (1-alfa)*dQdb(a(t),s(t),k) + alfa * tmp;
            qqb(:,:,:,t) = dQdb;
        end

		tmp = (dQde(:,s(t)) );
		dl(9) = dl(9) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-Q(a(t),s(t)))*alfa*(1-alfa);

		dl(10) = dl(10) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(11) = dl(11) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
	end

	Q(a(t),s(t)) = Q(a(t),s(t)) + alfa * (er - Q(a(t),s(t))); 
    qq(:,:,t)= Q;
	V(s(t))      = V(s(t))      + alfa * (er - V(s(t)     ));

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


