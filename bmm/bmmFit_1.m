function res = bmmFit_1(con,lb,ub)

rng default
nIter = 5; %number of initial search conditions


%param vector = [mu sd]
if ~exist('lb','var') || isempty(lb)
    lb = [0 0]';
end
if ~exist('ub','var') || isempty(ub)
    ub = [1 400]';
end
ubSet = [1 100]';


maxfeval = 10000;
options = optimset('TolX',1e-8,'maxfunevals',maxfeval, ...
    'maxiter',round(maxfeval/3),'GradObj','off','Display','notify', ...
    'DerivativeCheck','off','Algorithm','interior-point');

loglik = inf;
p = NaN;
for iIter = 1:nIter
    %fprintf(1,'iter %d of %d\n',iIter,nIter);
    x0 =  rand(numel(lb),1) .* (ubSet-lb) + lb;
    
    [x, ll] = fmincon(@(x) ll_bmmPDF_1(con,x), x0, ...
        [],[],[],[],lb,ub,[],options);
    
    if ll < loglik
        loglik = ll;
        p = x;
    end
end

res.nLL      = loglik; %negative log likelihood
res.w       = p(1);
res.c       = p(2);

