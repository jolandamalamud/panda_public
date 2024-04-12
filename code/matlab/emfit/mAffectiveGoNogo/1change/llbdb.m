function [l,dl,dsurr] = llbdb(x,D,mu,nui,doprior,options);
% 
% [l,dl,surrugatedata] = llb(x,D,mu,nui,doprior,options);
% 
% log likelihood (l) and gradient (dl) of only bias model allowing bias to
% change over sessions
% 
% Use this within emfit.m to tit RL type models to a group of subjects using
% EM. 

dodiff=nargout==2;
np = length(x);

% add Gaussian prior with mean mu and variance nui^-1 if doprior = 1 
[l,dl] = logGaussianPrior(x,mu,nui,doprior);

sess = vertcat(D.data.sess);
a = vertcat(D.data.a); 
r = vertcat(D.data.r); 
s = vertcat(D.data.s); 

bias = x(1);
delta = x(2);

rep = length(sess)/length(D.data);
rep = [1, rep+1, 2*rep+1];

if options.generatesurrogatedata==1
    pemp = ind_out_prob(s,a,r);
    a = zeros(size(a));
    dodiff=0;
end

for t=1:length(a)
    
    if ismember(t, rep)
       Q = init_Qvalue();
    end
    
    if isnan(a(t)); continue; end
    
    rho = sum(s(t)==[1 3]);

    q = Q(:,s(t)); 
    
    q(1) = q(1) + bias + sess(t) * delta;


    l0 = q - max(q);
    la = l0 - log(sum(exp(l0)));
    pa = exp(la); 

    if options.generatesurrogatedata==1
        [a(t),r(t)] = generatera(pa',s(t),pemp);
    end
    l = l + la(a(t));

    if dodiff
        dl(1) = dl(1) + (a(t)==1) - pa(1);
        dl(2) = dl(2) + sess(t) * ((a(t)==1) - pa(1));    
    end

    Q(a(t),s(t)) = Q(a(t),s(t));  

end
l  = -l; 
dl = -dl; 

if options.generatesurrogatedata==1
	dsurr.a = a; 
	dsurr.r = r; 
end


