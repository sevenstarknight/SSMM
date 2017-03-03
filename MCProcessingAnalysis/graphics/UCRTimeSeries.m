figure
ylim([-3.5, 6])
hold on

x = linspace(0, 1, length(structMCReduct(i).timeSeries(:,2)));

for i = 1:1:length(structMCReduct)
    
    
    if(structMCReduct(i).classType == 1)
        plot(x, structMCReduct(i).timeSeries(:,2), '-b')
    elseif(structMCReduct(i).classType == 2)
        plot(x, structMCReduct(i).timeSeries(:,2), ':g')
    elseif(structMCReduct(i).classType == 3)
        plot(x, structMCReduct(i).timeSeries(:,2), '--r')
    end
    
end

grid on
xlabel('Phase')
ylabel('(V- <V>)/\sigma_{V}')

legend('1', '2', '3')

hold off