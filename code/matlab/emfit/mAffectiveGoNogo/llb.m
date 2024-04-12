function [l,dl,dsurr] = llb(x,D,mu,nui,doprior,options)
% 
% [l,dl,surrugatedata] = llb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of random baseline model
% 
% Use this within emfit.m to tit RL type models to a group of subjects using
% EM. 

dodiff=nargout==2;
bias 		= x(1);						% constant bias 

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

a = D.a; 
r = D.r; 
s = D.s; 

Q=zeros(2,4); 

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
	a = zeros(size(a));
	dodiff=0;
end

for t=1:length(a)
    
    if a(t)==3; continue; end

	q = Q(:,s(t)); 
	q(1) = q(1) + bias;    % add Pavlovian effect 

	l0 = q - max(q);
	la = l0 - log(sum(exp(l0)));
	pa = exp(la); 

	if options.generatesurrogatedata==1
		[a(t),r(t)] = generatera(pa',s(t),pemp);
	end
	l = l + la(a(t));

    if dodiff
        dl(1) = dl(1) + (a(t)==1) - pa(1);
    end

	Q(a(t),s(t)) = Q(a(t),s(t));  

end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


