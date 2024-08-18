function [pdf, xc] = bmmPDF_1(w,c)
%clear;
%mu = [.5 .5];
%sd = .3;
%mix = .5; %wt assigned to first n-1 distributions

xe = 0:0.01:1;
dt = diff(xe(1:2));
xc = xe(1:end-1) + dt/2;

pdf1 = ( xc .^ ( c .* w(1) ) .* ( 1 - xc ) .^ ( c .* ( 1 - w(1) ) ) ) ./ ...
    beta( 1 + c .* w(1), 1 + c .* ( 1 - w(1) ) );


pdf = pdf1;
xc = xc';