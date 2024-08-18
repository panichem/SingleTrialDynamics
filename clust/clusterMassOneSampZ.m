function [pvals, clustMass, clustIdx, rawP] = clusterMassOneSampZ(x,varargin)

[nShuff, nTP] = size(x);

dat = x(1,:);
shuff = x(2:end,:);
z = ( dat - mean(shuff,1) ) ./ std(shuff,[],1);
mask = abs(z) > 1.96;
rawP = 1-normcdf(z);

[clustIdx, ~, clustMass] = getClust(mask,z);

nIter = nShuff-1;
nullMass = nan(nIter,1);
for iIter = 1:nIter
    
    idx = false(nShuff,1);
    idx(iIter+1) = true;

    tmpDat = x(idx,:);
    tmpShuff = x(~idx,:);
    z = ( tmpDat - mean(tmpShuff,1) ) ./ std(tmpShuff,[],1);
    mask = abs(z) > 1.96;

    [~, ~, tmpMass] = getClust(mask,z);

    if isempty(tmpMass)
        nullMass(iIter) = 0;
    else
        nullMass(iIter) = max(abs(tmpMass));
    end
end

pvals = sum(nullMass>=abs(clustMass)')./nIter;

