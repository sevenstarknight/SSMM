echo off
close all
clc
disp(' ')
disp('  This is a short demo of interval ECVA on spectroscopic data')
disp(' ')
disp('  We want to make an interval ECVA model of mid infrared data on 273 milk samples from four different animal groups')
disp('  The four groups are cow, goat, sheep and buffalo (data are provided by FOSS, www.foss.dk)')
disp(' ')
disp('  The spectroscopic range measured is from 932.91 (variable 1) to 2968.35 cm-1')
disp('  Due to intense water absorption the bands from 1580.55-1715.48 cm-1 and 1811.85-2698.50 cm-1 are removed')
disp(' ')
disp('  Load the data and plot the mean centered data')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
load DemoInfraredMilk
Xmncn = (X - ones(size(X,1),1)*mean(X)); % Mean center the data
figure(1)
% h=plot(wavenumber,Xmncn(1:100,:),'b',wavenumber,Xmncn(101:161,:),'g',wavenumber,Xmncn(162:233,:),'r',wavenumber,Xmncn(234:273,:),'c');
h=plot(1:266,Xmncn(1:100,:),'b',1:266,Xmncn(101:161,:),'g',1:266,Xmncn(162:233,:),'r',1:266,Xmncn(234:273,:),'c');
axis tight
title('Mean centered infrared spectra on milk samples')
xlabel('932.91 (variable 1) to 2968.35 cm-1')
ylabel('Absorbance')
legend(h([1 101 162 234]),'Cow    ','Goat   ','Sheep  ','Buffalo');
clear Xmncn
disp(' ')
disp('  Some spectal ranges seem relevant for discrimination but no range is capable of completely')
disp('  separation as judged from the inspection') 
disp(' ')
disp('  An interval ECVA model on the data with 12 intervals and a maximum of 15 PLS components in each is calculated')
disp('  using segmented systematic cross validation with ten segments')
disp(' ')
disp('  Use  iECVAm=ecva(''iecva'',X,classmarker,no_of_comp,prepro_method,intervals,xaxislabels,val_method,segments,prior,plots,minusglobal);')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
close
pause
echo on
iECVAm=ecva('iecva',X,classmarker,15,'none',12,wavenumber,'syst123',10);
echo off
disp(' ')
disp('  It is observed that interval 7 and 8 have 0 misclassifications while the global models has one error')
disp('  (it should be noted that if you use 17 PLS components in the global model the number of misclassifications is zero)')
disp(' ')
disp('  We would now like to see the same plot with variable number instead of interval numbers')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
ecva('iplot',iECVAm,'varlabel');
echo off
disp(' ')
disp('  Finally we will inspect the canonical variates plots for each interval; 1:Cow, 2:Goat, 3:Sheep, 4:Buffalo')
disp('  Please note that you can move legends that overlap the data points by just dragging them away from the plot')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
ecva('icanonvar',iECVAm,1,2,[]);
echo off
disp(' ')
disp('  END OF DEMO')
