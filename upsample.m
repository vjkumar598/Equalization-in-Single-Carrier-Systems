% upsample: i.e, increase the sampling rate of the given vector x by a
% factor m (positive integer)
function xup = upsample(x,m)
    
    lus = (length(x)- 1)*m + 1; % length of upsampling vector
    xup = zeros(1,lus);
    
    xup(1:m:lus) = x;
end