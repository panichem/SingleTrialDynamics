function ll = ll_bmmPDF(con,x)
%output negative log likelihood           

%parse params
w = x(1:2);
c = x(3); 
mix = x(4);

[P, xc] = bmmPDF(w,c,mix);
[~, i] = min(abs(xc - con'),[],1); 
p = P(i);
ll = sum(-log(p));
if ~isreal(ll); ll = inf; end