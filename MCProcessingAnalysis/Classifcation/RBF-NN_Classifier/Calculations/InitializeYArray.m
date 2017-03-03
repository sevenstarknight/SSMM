function [y] = InitializeYArray(structLRC,grpTraining)

uniqueGrps = unique(grpTraining);

intUniqueGroups = structLRC.intUniqueGroups;
intLengthTraining = structLRC.intLengthTraining;

y = zeros(intLengthTraining,intUniqueGroups);

for i = 1:1:intLengthTraining
    
    for j = 1:1:intUniqueGroups
        if( strcmp(grpTraining(i), uniqueGrps(j)) ) 
            y(i,j) = 1;
        else
            y(i,j) = 0;
        end
    end
    
end

end