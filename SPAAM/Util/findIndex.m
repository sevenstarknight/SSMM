function[index] = findIndex(sample, states)

    if(sample >= max(states))
        index = length(states);            
    else
        %%lower bound
        index = find(sample < states, 1) - 1;        
        if(index == 0) 
            index = 1;
        end
    end 

end