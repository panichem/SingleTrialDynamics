function ll = ll_bmmPDF_1(con,x)
%output negative log likelihood           

%parse params
w = x(1);
c = x(2); 

[P, xc] = bmmPDF_1(w,c);
[~, i] = min(abs(xc - con'),[],1); 
p = P(i);
ll = sum(-log(p));
if ~isreal(ll); ll = inf; end