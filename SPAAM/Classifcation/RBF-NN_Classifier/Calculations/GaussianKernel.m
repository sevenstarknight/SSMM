function [fDistance] = GaussianKernel(fltPositionVector, nodes, spread)

fDistance = zeros(1,length(nodes));
intDimensions = length(fltPositionVector);

for i = 1:1:length(nodes)
    centerEstimate = nodes{i};
    distance = fltPositionVector - centerEstimate;
    
    if(iscell(spread))
        cov = spread{i};
        detCov = 1/sqrt(det(cov));
        
        fltTDistance = ((distance/cov)*distance');
        
        fDistance(i) = ((2*pi)^(-intDimensions/2))*detCov*exp(-0.5*fltTDistance);
    else
        fDistance(i) = exp(-(norm(distance))/spread);
    end
end

end