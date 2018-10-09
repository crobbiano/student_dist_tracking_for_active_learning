function [ newdist1, newdist2, bins ] = makeSimilarDists( dist1, bins1, dist2, bins2, binwidth)
    %makeSimilarDists Makes dists1 and dists2 span the same bin edges
    %   Makes new bins from min(min(dist1), min(dist2)) to 
    %   max(max(dist1), max(dist2)) for newdist1 and newdist2 and places 
    %   the values from dist1,2 into the newdist1,2 at the correct bin
    %   This is so we can run distance metrics on them to compare the dists
    
    binmax = max(max(bins1), max(bins2));
    binmin = min(min(bins1), min(bins2));
    
    bins = binmin:binwidth:binmax;
    
    newdist1 = zeros(1, numel(bins));
    for i=1:numel(dist1)  
        [val, idx] = find(abs(bins - bins1(i)) < 10^-3);
        newdist1(idx) = dist1(i);
    end
    
    newdist2 = zeros(1, numel(bins));
    for i=1:numel(dist2)
        [val, idx] = find(abs(bins - bins2(i)) < 10^-3);
        newdist2(idx) = dist2(i);
    end
end

