function [track, l,dl,dsurr] = ll2aepxbcon(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llbaepxb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of simple RW model with constant bias
% towards one action, irreducible noise and positive Pavlovian bias parameter
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
alfa 		= 1./(1+exp(-x(1:2)));		% learning rate
epsilon 	= exp(x(3));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
g       	= 1/(1+exp(-x(4)));		% irreducible noise
bias 		= x(5);						% constant bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

VP = zeros(4,1);
VI = zeros(2,4);
MQ = alfa(1) + zeros(2,4);
MV = alfa(2) + zeros(4,1);
L = zeros(4,1);

dQdb = zeros(2,4);
dQde = zeros(2,4);
dVdb = zeros(4,1);
dVde = zeros(4,1);
dwdb = 0; dLdb = 0;

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
    
%     if t == 1; L = log(w0) - log(1-w0); end

    w = 1./(1+exp(-L(s(t))));
    q = [epsilon * w * VP(s(t)) + bias;0] + (1 - w) * VI(:,s(t));
    
    l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
    end
    
	l = l + log(pg(a(t)));

% 	er = beta * r(t);
% 
% 	if dodiff
% 		tmp = (1 - w) * dQdb(:,s(t)) + (1 - dwdb) * VI(:,s(t)) + [epsilon*w*dVdb(s(t)) + epsilon*dwdb*VP(s(t));0];
% 		dl(1) = dl(1) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
%         dwdb = 1./(1+exp(-dLdb));
%         dLdb = dLdb + r(t) * (safelog(dVdb(s(t))) - safelog(dQdb(a(t),s(t)))) + ...
%             (1 - r(t)) * (safelog(1 - dVdb(s(t))) - safelog(1 - dQdb(a(t),s(t))));
%         dQdb(a(t),s(t)) = (1-alfa)*dQdb(a(t),s(t)) + alfa*er;
% 		dVdb(s(t))      = (1-alfa)*dVdb(s(t))      + alfa*er;
% 
% 		tmp = (1 - w) * dQde(:,s(t)) + [epsilon*w*dVde(s(t));0];
% 		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
% 		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-VI(a(t),s(t)))*alfa*(1-alfa);
% 		dVde(s(t))      = (1-alfa)*dVde(s(t))      + (er-VP(s(t)))*alfa*(1-alfa);
% 
% 		dl(3) = dl(3) + g*(p0(a(t))*epsilon*w*VP(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));
% 
% 		dl(4) = dl(4) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));
% 
% 		tmp = [1;0];
% 		dl(5) = dl(5) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
%     end   
    
%      if sum(s(t) == [1,3])
%         if r(t) == 1
%             update = r(t) * (log(VP(s(t))) - log(VI(a(t),s(t))));
%         else
%             update = (1 - r(t)) * (log(1 - VP(s(t))) - log(1 - VI(a(t),s(t))));
%         end
%     else
%         if r(t) == -1
%             update = - r(t) * (log(-VP(s(t))) - log(-VI(a(t),s(t))));
%         else
%             update = (1 + r(t)) * (log(1 + VP(s(t))) - log(1 + VI(a(t),s(t))));
%         end
%     end

    update = log(VP(s(t))/VI(a(t),s(t)));

    if ~isnan(update) && ~isinf(update); L(s(t)) = L(s(t)) + update; end
    
	VI(a(t),s(t))   = VI(a(t),s(t)) + alfa(1) * (r(t) - VI(a(t),s(t)));
	VP(s(t))        = VP(s(t))      + alfa(2) * (r(t) - VP(s(t)));
    
    track.Q(:,:,t) = VI; track.V(:,t) = VP; track.L(:,t) = L; track.w(t) = w; track.u(t) = update;

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

