

clc
clear

stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab\data\dataset';

cd(stringDirectory);

listing = dir(stringDirectory);

listing = listing(3:end);

errorArray1 = zeros(1,length(listing));
for i = 1:1:length(listing)
    strLabel = listing(i).name;
    disp(strcat(num2str(i), ': ', strLabel))
    errorArray(i) = SPAAM(stringDirectory, strLabel);
end

stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab';
cd(stringDirectory);

save('errorForSSMM_1', 'errorArray');


stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab\data\dataset2';

cd(stringDirectory);

listing2 = dir(stringDirectory);
listing2 = listing2(3:end);

errorArray2 = zeros(1,length(listing));
for i = 1:1:length(listing2)
    strLabel = listing2(i).name;
    disp(strcat(num2str(i), ': ', strLabel))
    errorArray2(i) = SPAAM(stringDirectory, strLabel);
end

stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab';
cd(stringDirectory);

save('errorForSSMM_2', 'errorArray2');