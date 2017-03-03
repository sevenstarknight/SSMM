PLEASE READ THESE FEW LINES BEFORE STARTING:

MATLAB should be installed. No other toolboxes are needed.

The main file is the ecva.m file from which you can perform several different tasks. Just type ecva followed by enter to get help.

ecva('model'): You can see the number of misclassifications by typing: ECVAm.MisClassVal (assuming that you call the output for ECVAm).

Get the class predictions by typing: ECVAm.ClassPredOptComp.

ecva('cv1d')/ecva('cv2d') and ecva('cw1d')/ecva('cw2d') plot the results.

ecva('pred') is for prediction based on results from the ecva('model').

If you work on e.g. chemical data (different units) you can use 'scale' as the prepro_method instead of 'none'.
NOTE: the ECVA method does not use mean centering or autoscaling as e.g. PCA and PLS methods, but 'none' or 'scaling'.

Have a look at ecvademo1.m and ecvademo2.m to get inspired.

ecva('iecva') calculates interval ECVA models and ecva('iplot')/ecva('icanonvar') plots the results.

More information can be found in:
Nørgaard L, Bro R, Westad F, Engelsen SB. A modification of Canonical Variates Analysis to handle highly collinear multivariate data, Journal of Chemometrics, 20:425-435, 2006.

In the paper it is stated that quadratic LDA is also possible (not a true quadratic ECVA, only the LDA part is quadratic). This is not implemented in the present version of the toolbox.

=============== END ===============
