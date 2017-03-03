echo off
close all
clc
disp(' ')
disp('  This is a short demo of ECVA')
disp(' ')
disp('  We want to make an ECVA model of fluorescence spectra')
disp('  on 83 samples from six different factory sites')
disp(' ')
disp('  Load the data')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
load DemoFluoData
echo off
disp(' ')
disp('  Let''s plot the concatenated (and 11-data points smoothed) spectra')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
plot(1:1023,X(1:13,:),'k',1:1023,X(14:27,:),'b',1:1023,X(28:42,:),'r',1:1023,X(43:56,:),'c',1:1023,X(57:71,:),'m',1:1023,X(72:83,:),'g')
axis tight
echo off
disp(' ')
disp('  The model is run with 20 components')
disp('  and 5 segments using syst123 cross validation')
disp(' ')
disp('  Use ECVAm = ecva(''model'',X,classmarker,no_of_lv,prepro_method,cv_method,segments,prior,plots);')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
ECVAm = ecva('model',X,classmarker,20,'none','syst123',5);
echo off
disp(' ')
disp('  The two plots show the canonical variates (CV) and the corresponding weights for the data')
disp('  An almost perfect separation of the six groups is observed in the CV plot')
disp('  showing the two first of the five CVs. Using all five CVs only one sample')
disp('  is misclassified in the cross validated model.')
disp(' ')
disp('  END OF DEMO')
