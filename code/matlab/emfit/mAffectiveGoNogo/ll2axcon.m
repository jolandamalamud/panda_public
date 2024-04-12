function [l,dl,dsurr] = ll2axcon(x,D,mu,nui,doprior,options);
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
alfa        = 1./(1+exp(-x(1:2)));
g       	= 1/(1+exp(-x(3)));		% irreducible noise

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

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
    
%     if t == 1; L = log(w0) - log(1-w0); end

    w = 1./(1+exp(-L(s(t))));
    q = w * VP(s(t)) + (1 - w) * VI(:,s(t));
    
    l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
    end
    
	l = l + log(pg(a(t)));

% 	if dodiff
%         tmp = (1 - w) * dQde(:,s(t)) +  (1 - dwe) * VI(:,s(t));
% 		dl(1) = dl(1) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
%         dwe = 1./(1+exp(-dLe));
%         dLe = 
%         if ~isnan(update) && ~isinf(update); dLe = dLe + update; end
% 		dQde(a(t),s(t)) = (1-alfa(1))*dQde(a(t),s(t)) + (r(t)-VI(a(t),s(t)))*alfa(1)*(1-alfa(1));
%         tmp = [epsilon*w*dVde(s(t));0];
% 		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
% 		dVde(     s(t)) = (1-alfa(2))*dVde(     s(t)) + (r(t)-V(s(t)))*alfa(2)*(1-alfa(2));
% 
% 		tmp = (1 - w) * dQde(:,s(t)) + [epsilon*w*dVde(s(t));0];
% 		dl(2) = dl(2) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
% 		dQde(a(t),s(t)) = (1-alfa)*dQde(a(t),s(t)) + (er-VI(a(t),s(t)))*alfa*(1-alfa);
% 		dVde(s(t))      = (1-alfa)*dVde(s(t))      + (er-VP(s(t)))*alfa*(1-alfa);
% 
% 		tmp = [1;0];
% 		dl(5) = dl(5) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
%     end   
    
    if sum(s(t) == [1,3])
        if r(t) == 1
            update = r(t) * (log(VP(s(t))) - log(VI(a(t),s(t))));
        else
            update = (1 - r(t)) * (log(1 - VP(s(t))) - log(1 - VI(a(t),s(t))));
        end
    else
        if r(t) == -1
            update = - r(t) * (log(-VP(s(t))) - log(-VI(a(t),s(t))));
        else
            update = (1 + r(t)) * (log(1 + VP(s(t))) - log(1 + VI(a(t),s(t))));
        end
    end

%     update = log(VP(s(t))/VI(a(t),s(t)));

    if ~isnan(update) && ~isinf(update); L(s(t)) = L(s(t)) + update; end
    
%     MQ(a(t),s(t)) = MQ(a(t),s(t)) + 1;
%     MV(s(t)) =  MV(s(t)) + 1;
	VI(a(t),s(t))   = VI(a(t),s(t)) + alfa(1) * (r(t) - VI(a(t),s(t)));
	VP(s(t))        = VP(s(t))      + alfa(2) * (r(t) - VP(s(t)));
    
%     track.Q(:,:,t) = VI; track.V(:,t) = VP; track.L(:,t) = L; track.w(t) = w;

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

