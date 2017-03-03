function [structProb] = DetermineDistanceEstimates(fltPattern, structParaEst)

structProb = struct([]);

for i = 1:1:length(structParaEst)
    tmpStruct = structParaEst(i);
    
    mu = tmpStruct.mean;
    invCov = tmpStruct.invCov;
    constant = tmpStruct.constant;

    structProb(i).prob = constant - 0.5*(fltPattern-mu)*invCov*(fltPattern-mu)';
    
    structProb(i).type = tmpStruct.type;
end
    
end

