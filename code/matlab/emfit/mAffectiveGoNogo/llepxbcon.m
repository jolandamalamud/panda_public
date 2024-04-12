function [track,l,dl,dsurr] = llepxbcon(x,D,mu,nui,doprior,options);
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
epsilon 	= exp(x(1));				% 'pavlovian' parameter. Weigth of Vcue into Qgo
g       	= 1/(1+exp(-x(2)));		% irreducible noise
bias 		= x(3);						% constant bias 
alfa        = 1/(1+exp(-x(4:5)));

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

VP = zeros(4,1);
VI = zeros(2,4);
MQ = zeros(2,4);
MV = zeros(4,1);
L = zeros(4,1);

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
    
%     if t == 1; L = log(w0) - log(1-w0); end

    w = 1./(1+exp(-L(s(t))));
%     q = [epsilon * w * VP(s(t)) + bias;0] + (1 - w) * VI(:,s(t));
    q = w * VP(s(t)) + (1 - w) * VI(:,s(t));
    
    l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	p0 = exp(la); 
	pg = g*p0 + (1-g)/2;

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pg',s(t),pemp);
    end
    
	l = l + log(pg(a(t)));

	if dodiff
		dl(1) = dl(1) + g*(p0(a(t))*epsilon*w*VP(s(t)) * ((a(t)==1)-p0(1))) / pg(a(t));

		dl(2) = dl(2) + g*(1-g)*(p0(a(t))-1/2)/pg(a(t));

		tmp = [1;0];
		dl(3) = dl(3) + g*(p0(a(t)) * (tmp(a(t)) - p0'*tmp)) / pg(a(t));
    end   
    
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
    disp(update)

%     update = log(VP(s(t))/VI(a(t),s(t)));

    if ~isnan(update) && ~isinf(update); L(s(t)) = L(s(t)) + update; end
    
    MQ(a(t),s(t)) = MQ(a(t),s(t)) + 1;
    MV(s(t)) =  MV(s(t)) + 1;
	VI(a(t),s(t))   = VI(a(t),s(t)) + alfa * (r(t) - VI(a(t),s(t)))/MQ(a(t),s(t));
	VP(s(t))        = VP(s(t))      + alfa * (r(t) - VP(s(t)))/MV(s(t));
    
    track.Q(:,:,t) = VI; track.V(:,t) = VP; track.L(:,t) = L; track.w(t) = w;

end
l  = -l;
dl = -dl;

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end

