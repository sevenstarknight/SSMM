echo off
close all
clc
disp(' ')
disp('  This is a short demo of ECVA on chemical data')
disp(' ')
disp('  We want to make an ECVA model of quality measurements')
disp('  on 41 samples from three different factory sites')
disp(' ')
disp('  Load the data and plot the autoscaled chemical data in a line plot')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
load DemoChemicalData
Xauto = (X - ones(size(X,1),1)*mean(X))./(ones(size(X,1),1)*std(X)); % Autoscales the data
figure(1)
h=plot(1:10,Xauto(1:15,:),'b',1:10,Xauto(16:29,:),'g',1:10,Xauto(30:41,:),'r');
axis tight
legend(h([1 16 30]),'Factory 1','Factory 2','Factory 3'); clear h
clear Xauto
disp(' ')
disp('  The ''green'' factory is clearly different from the ''red'' and ''blue''') 
disp('  at the first variables while ''red'' and ''blue'' are strongly overlapping')
disp(' ')
disp('  An ECVA model on scaled data with 10 components')
disp('  and full cross validation is now calculated')
disp(' ')
disp('  Use ECVAm = ecva(''model'',X,classmarker,no_of_lv,prepro_method,cv_method,segments,prior,plots);')
disp(' ')
disp('  Press any key to continue (we suppress plots)')
disp(' ')
pause
echo on
ECVAm = ecva('model',X,classmarker,10,'scale','syst123',10,[],0);
echo off
disp(' ')
disp('  Now plot the extended canonical variates (ECVs) [Press any key to continue]')
disp(' ')
pause
echo on
figure(2)
ecva('cv2d',ECVAm,1,2)
echo off
disp(' ')
disp('  A clear separation is observed in the CV plot')
disp(' ')
disp('  The weights are plotted in a line plot. If we want to see a 2D plot we use:')
disp('  ecva(''cw2d'',w_a,w_b,VarNames,no_of_comp)')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
figure(3)
ecva('cw2d',ECVAm,1,2,VarNames)
echo off
disp(' ')
disp('  For comparison we will make a score plot from an autoscaled PCA model')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
figure(4)
subplot(121)
ecva('cv2d',PCAmodel,1,2) % Plots PCA scores
xlabel('PC#1')
ylabel('PC#2')
title('PCA')
subplot(122)
ecva('cv2d',ECVAm,1,2)
title('ECVA')
echo off
clear PCAmodel
disp(' ')
disp('  It is seen that the ECVA method produces a more discriminative ''score'' plot')
disp(' ')
disp('  END OF DEMO')
