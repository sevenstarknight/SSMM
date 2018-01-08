x = linspace(0, 1, length(structMCReduct(1).timeSeries(:,2)));

figure()
hold on
for indx = idxTrainCross
    
    if(strcmp(grpSource{indx},'1'))
        plot(x, structMCReduct(indx).timeSeries(:,2), '.r')
    elseif(strcmp(grpSource{indx},'2'))
        plot(x, structMCReduct(indx).timeSeries(:,2), ':b')
    elseif(strcmp(grpSource{indx},'3'))
        plot(x, structMCReduct(indx).timeSeries(:,2), '-g')
    end
       
end

xlabel('Phase')
ylabel('Amplitude')

hold off

%% ====================================

figure()
hold on
count = 1;
for indx = idxTrainCross
    
    if(intLabel(count) == 1)
        plot(x, structMCReduct(indx).timeSeries(:,2), '-r')
    elseif(intLabel(count) == 2)
        plot(x, structMCReduct(indx).timeSeries(:,2), ':b')
    elseif(intLabel(count) == 3)
        plot(x, structMCReduct(indx).timeSeries(:,2), '-g')
    elseif(intLabel(count) == 4)
        plot(x, structMCReduct(indx).timeSeries(:,2), ':k')
    elseif(intLabel(count) == 5)
        plot(x, structMCReduct(indx).timeSeries(:,2), '-m')
    end
       
    count = count + 1;
end

xlabel('Phase')
ylabel('Amplitude')

%% ====================================

figure()
hold on
count = 1;
for indx = idxTrainCross
    
    if(intLabel(count) == 4)
        plot(x, structMCReduct(indx).timeSeries(:,2), '-r')
    end
       
    count = count + 1;
end

xlabel('Phase')
ylabel('Amplitude')