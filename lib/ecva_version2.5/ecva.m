function varargout = ecva(command,varargin)

%  The ECVA Toolbox
%     Version 2.5 - November 2010
%     Copyright(c) Nørgaard & Bro, Chemometrics Group, University of Copenhagen, Denmark
%
%  For help on how to use ECVA for different purposes, type
%
%  ECVA MODELS
%     ecva('model'): Calculates an ECVA model
%     ecva('pred'):  Used for classification of new samples
%     ecva('iecva'): Calculates an interval ECVA model
%
%  ECVA PLOTS
%     ecva('cv1d'):  Bar plot of canonical variate(s)
%     ecva('cv2d'):  2D-plot of canonical variates
%     ecva('cw1d'):  1D-plot canonical weights
%     ecva('cw2d'):  2D-plot canonical weights
%     ecva('iplot'): Plot of iECVA results
%     ecva('icanonvar'): Plot iECVA canonical variates
%
%  MISC FILES
%     ecva('intervals'): Table with information on the intervals in an iECVA model
%
%  ECVA DEMOS
%     ecva('demo1'): ECVA demo on near infrared spectroscopic data
%     ecva('demo2'): ECVA demo on chemical data
%     ecva('demo3'): interval ECVA demo on infrared spectroscopic data
%
%  ECVA TEST FOR OVERFIT
%     ecva('randtest'): Tests overfit for a given data dimension

if nargin==0 % Help
    disp(' ')
    disp('   The ECVA Toolbox')
    disp('      Version 2.5 - November 2010')
    disp('      Copyright(c) Nørgaard & Bro, Chemometrics Group, University of Copenhagen, Denmark')
    disp(' ')
    disp('   For help on how to use ECVA for different purposes, type')
    disp(' ')
    disp('   ECVA MODELS')
    disp('      ecva(''model''): Calculates an ECVA model')
    disp('      ecva(''pred''):  Used for classification of new samples')
    disp('      ecva(''iecva''): Calculates an interval ECVA model')
    disp(' ')
    disp('   ECVA PLOTS')
    disp('      ecva(''cv1d''):  Bar plot of canonical variate(s)')
    disp('      ecva(''cv2d''):  2D-plot of canonical variates')
    disp('      ecva(''cw1d''):  1D-plot canonical weights')
    disp('      ecva(''cw2d''):  2D-plot canonical weights')
    disp('      ecva(''iplot''): Plot of iECVA results')
    disp('      ecva(''icanonvar''): Plot iECVA canonical variates')
    disp(' ')
    disp('   MISC FILES')
    disp('      ecva(''intervals''): Table with information on the intervals in an iECVA model')
    disp(' ')
    disp('   ECVA DEMOS')
    disp('      ecva(''demo1''): ECVA demo on near infrared spectroscopic data')
    disp('      ecva(''demo2''): ECVA demo on chemical data')
    disp('      ecva(''demo3''): interval ECVA demo on infrared spectroscopic data')
    disp(' ')
    disp('   ECVA TEST FOR OVERFIT')
    disp('      ecva(''randtest''): Tests overfit for a given data dimension')
    disp(' ')

elseif nargin == 1 % Help
    switch lower(command)
        case {'model'}
            help ecva>a_ecvamodel
        case {'pred'}
            help ecva>b_ecvapred
        case {'cv1d'}
            help ecva>c1_cv1d
        case {'cv2d'}
            help ecva>c2_cv2d
        case {'cw1d'}
            help ecva>d_cw1d
        case {'cw2d'}
            help ecva>e_cw2d
        case {'iecva'}
            help ecva>f_iecva
        case {'iplot'}
            help ecva>g_iplot
        case {'icanonvar'}
            help ecva>h_icanonvar
        case {'intervals'}
            help ecva>i_intervals
        case {'demo1'}
            ecvademo1; % Separate m-file
        case {'demo2'}
            ecvademo2; % Separate m-file
        case {'demo3'}
            ecvademo3; % Separate m-file
        case {'randtest'}
            help ecva>k_randtest
        otherwise
            disp('Not valid command')
    end
else % Modelling
    switch lower(command)
        case {'model'}
            varargout{:} = a_ecvamodel(varargin{:});
        case {'pred'}
            varargout{:} = b_ecvapred(varargin{:});
        case {'cv1d'}
            c1_cv1d(varargin{:});
        case {'cv2d'}
            c2_cv2d(varargin{:});
        case {'cw1d'}
            d_cw1d(varargin{:});
        case {'cw2d'}
            e_cw2d(varargin{:});
        case {'iecva'}
            varargout{:} = f_iecva(varargin{:});
        case {'iplot'}
            g_iplot(varargin{:});
        case {'icanonvar'}
            h_icanonvar(varargin{:});
        case {'intervals'}
            i_intervals(varargin{:});
        case {'demo1'}
            ecvademo1;
        case {'demo2'}
            ecvademo2;
        case {'demo3'}
            ecvademo3;
        case {'randtest'}
            varargout{:} = k_randtest(varargin{:});
        otherwise
            disp('Not valid command')
    end
end

function varargout = a_ecvamodel(varargin)
%
%  Extented Canonical Variates Analysis (ECVA) with cross validation
%
%  Input:
%    X: a double array containing the independent variables
%    classmarker (column vector): contains the group information (you can use either numbers or strings).
%    no_of_comp: the number of PLS components in the inner CVA relation
%              NOTE: choose a large number since this PLS model includes a variables x variables matrix
%              The optimal number will be found by cross validation (if this option is selected)
%    prepro_method: 'none' (e.g. NIR spectra) or 'scale' (1/Sdev; e.g. chemical data)
%    cv_method: 'full' (leave-one-out), 'syst111' (cont. blocks), 'syst123' (Venetian blinds), 'random', 'manual' or 'none'
%               if 'none' is selected, insert [] in segments (or only use five inputs)
%    segments: segments in cross validation
%              Using the number of samples in X corresponds to full cv ('full' overrules segments)
%              if cv_methods is 'manual' use script_manseg to form the segments,
%              then insert the cell array formed (called manualseg) instead of the number of segments
%    prior: [optional], if set to [] equal group probability will be used
%    plots: [optional], if set to 0 no plots & text will be output, otherwise set to 1
%
%  Output:
%    ECVAm: a structured array containing all model and validation information
%           ECVAm.MisClassVal contains the number of misclassifications
%             (validated) as a function of the number of PLS components
%           ECVAm.ClassPredOptComp to see the class predictions for the
%             optimal number of PLS components
%
%  ECVAm = ecva('model',X,classmarker,no_of_comp,prepro_method,cv_method,segments,prior,plots);
%  Examples: 
%  ECVAm = ecva('model',X,classmarker,10,'none','syst123',15,[],0);
%    or
%  ECVAm = ecva('model',X,classmarker,10,'none','syst123',15,[0.25 0.25 0.50]);
%    or
%  ECVAm = ecva('model',X,classmarker,10,'scale','syst123',15);

X=varargin{1};classmarker=varargin{2};no_of_comp=varargin{3};prepro_method=varargin{4};cv_method=varargin{5};
plots=1; % Plots & text as output [default]
switch nargin
    case 6
        segments=varargin{6};
    case 7
        segments=varargin{6};
        prior=varargin{7};
    case 8
        segments=varargin{6};
        prior=varargin{7};
        plots=varargin{8};        
end

if ~ismember(prepro_method,{'none', 'scale'})
    disp(' ')
    disp('   Use ''none'' or ''scale'' as input for prepro_method')
    disp('   prepro_method is set to ''none''')
    disp(' ')
    prepro_method = 'none';
end

ECVAm=struct('MisClassVal',[],'MisClassCal',[],'OptComp',[],'ClassPredOptComp',[],'ClassProbOptComp',[],...
                    'CanonicalWeights',[],'CanonicalVariates',[],'ValidationMethod',[],'PreproMethod',[],'Prior',[],'ClassmarkerOriginal',[],'ClassmarkerNumeric',[],'Detail',[],'Compression',[]);
ECVAm.PreproMethod=prepro_method;
ECVAm.ValidationMethod=cv_method;

ECVAm.Detail.mx = mean(X); % To be used for predictions
ECVAm.Detail.stdx = std(X); % To be used for predictions
ECVAm.Compression = 'Yes';
        
% Data compression
[NoofSamples,NoofVariables] = size(X);
if NoofSamples < NoofVariables
    if strcmpi(prepro_method,'scale')
        % Note: in next version this should be integrated properly in the cross validation part        
        X = X./(ones(NoofSamples,1)*std(X));
    end
    [U,S,V] = svd(X,'econ');
    X = U*S;
    ECVAm.Compression = 'Yes';
end

X=sub_makestruct(X,classmarker);

ECVAm.ClassmarkerOriginal=X.ClassmarkerOriginal;
ECVAm.ClassmarkerNumeric=X.ClassmarkerNumeric;

if X.noofgroups > size(X.raw,2)
    min_no_of_comp=size(X.raw,2);
    max_no_of_comp=min_no_of_comp;
    if no_of_comp~=max_no_of_comp
        s=sprintf('Number of PLS components is locked to %g',max_no_of_comp);
        disp(' ')
        disp(s)
    end
else
    min_no_of_comp=X.noofgroups; % Or only min_no_of_comp=X.noofgroups
    max_no_of_comp=size(X.raw,2);
    if no_of_comp < min_no_of_comp
        s=sprintf('Maximum number of PLS components is set to %g',min_no_of_comp);
        disp(' ')
        disp(s)
        max_no_of_comp = X.noofgroups;
    elseif no_of_comp > max_no_of_comp
        s=sprintf('Maximum number of PLS components is set to %g',max_no_of_comp);
        disp(' ')
        disp(s)
    else
        max_no_of_comp = no_of_comp;
    end
end

if plots ~= 0
    messagetext=sprintf('There are %g unique classes',X.noofgroups);
    disp(' ')
    disp(messagetext)
    disp(' ')
end

if nargin <= 6 || isempty(prior)
    prior=ones(1,X.noofgroups-1)/X.noofgroups;
    prior=[prior 1-sum(prior)];
end

ECVAm.Prior=prior;

switch cv_method
    case {'full', 'syst111', 'syst123', 'random', 'manual'}
        if nargin >= 6 && iscell(segments)
          manualseg=segments; % A cell array
          segments=max(size(segments)); % Now a scalar
        end

        n=size(X.raw,1);
        no_sampl=fix(n/segments);
        left_over_samples=mod(n,segments);

        if strcmpi(cv_method,'full')
           cv_method='syst123';
           segments=n;
        end

        if strcmpi(cv_method,'random')
           ix=randperm(n);
        end

        count=1;
        for i=1:segments
            if plots ~= 0
                s=sprintf('Cross validation segment number %g',i); disp(s)
            end
            if strcmpi(cv_method,'syst111')
                  if left_over_samples == 0
                       count=count;
                       p_cvs=((i-1)*no_sampl+1+(count-1):i*no_sampl+(count-1))';
                  else   
                       p_cvs=((i-1)*no_sampl+1+(count-1):i*no_sampl+count)';
                       count=count+1;
                       left_over_samples=left_over_samples-1;
                  end
             elseif strcmpi(cv_method,'syst123')
                  p_cvs=(i:segments:n)';
             elseif strcmpi(cv_method,'random')
                  p_cvs=(i:segments:n)';
                  p_cvs=ix(p_cvs)';
             elseif strcmpi(cv_method,'manual')
                if max(size(manualseg)) ~= segments
                   disp('The number of segments does not correspond to the segments in manualseg')
                   break
                end
                Nsamples=0;
                for j=1:segments
                   Nsamples=Nsamples + max(size(manualseg{j}));
                end
                if n ~= Nsamples
                   disp('The number of samples in X does not correspond to the number of samples in manualseg')
                   break
                end
                p_cvs=manualseg{i};
             end
             tot=(1:n)';
             tot(p_cvs)=[];
             m_cvs=tot;
             ECVAm.Detail.CrossValSegments{i}=p_cvs;
             Xseg=X.raw(m_cvs,:);
             Xpseg=X.raw(p_cvs,:);
             if strcmpi(prepro_method,'scale') && (NoofSamples >= NoofVariables)
                   stdx=std(Xseg);     
                   Xseg=Xseg./(ones(length(m_cvs),1)*stdx);
                   Xpseg=Xpseg./(ones(length(p_cvs),1)*stdx);
             end
             ECVAseg=sub_ecva(Xseg,X.ClassmarkerNumeric(m_cvs),min_no_of_comp,max_no_of_comp);
             % Mean center prediction objects before projection onto W
             [~,mx]=sub_mc(Xseg);
             Xpseg=Xpseg - ones(size(Xpseg,1),1)*mx;
             for j=min_no_of_comp:max_no_of_comp
                  PredScores{j}=Xpseg*ECVAseg.CVA_weights{j};
             end
             for j=min_no_of_comp:max_no_of_comp
                 [ClassPred{j}(p_cvs,:),ProbPred{j}(p_cvs,:)]=sub_lda(PredScores{j},ECVAseg.ProjScores{j},X.ClassmarkerNumeric(m_cvs),prior);
             end
        end
    case {'none'} % No validation
         Xnoval=X.raw;
         if strcmpi(prepro_method,'scale') && (NoofSamples >= NoofVariables)
            stdx=std(Xnoval);     
            Xnoval=Xnoval./(ones(size(X.raw,1),1)*stdx);
         end
         ECVAseg=sub_ecva(Xnoval,X.ClassmarkerNumeric,min_no_of_comp,max_no_of_comp);
         % Mean center calibration objects before projection onto W
         [Xnoval,~]=sub_mc(Xnoval);
         for j=min_no_of_comp:max_no_of_comp
            PredScores{j}=Xnoval*ECVAseg.CVA_weights{j}; % Autoprediction, i.e. fit
         end
         for j=min_no_of_comp:max_no_of_comp % To handle few variables
            [ClassPred{j},ProbPred{j}]=sub_lda(PredScores{j},ECVAseg.ProjScores{j},X.ClassmarkerNumeric,prior);
         end
    otherwise
         disp('Not accepted validation method')
end

ECVAm.Detail.ClassPred=ClassPred;

ECVAm.MisClassVal = NaN*zeros(1,max_no_of_comp);
for j=min_no_of_comp:max_no_of_comp
    ECVAm.MisClassVal(j)=length(X.ClassmarkerNumeric)-size(find((ClassPred{j} - X.ClassmarkerNumeric)==0),1);
end

[~,ECVAm.OptComp]=min(ECVAm.MisClassVal);
% From 'help min': If the values along the first non-singleton dimension contain more
% than one minimal element, the index of the first one is returned.

ClassmarkerOriginalUnique = unique(ECVAm.ClassmarkerOriginal);
ClassmarkerNumericUnique = unique(ECVAm.ClassmarkerNumeric);
for i = 1:length(ClassPred{ECVAm.OptComp})
    ixtemp = ClassPred{ECVAm.OptComp}(i) == ClassmarkerNumericUnique;

    try
        ECVAm.ClassPredOptComp(i) = ClassmarkerOriginalUnique(ixtemp); % Does this also work if original classmarker is char?
    catch me       
       ixtemp; 
    end
end
ECVAm.ClassPredOptComp=ECVAm.ClassPredOptComp';
ECVAm.ClassProbOptComp=ProbPred{ECVAm.OptComp};
if ~isnumeric(ECVAm.ClassmarkerOriginal)
    ECVAm.ClassPredOptComp = char(ECVAm.ClassPredOptComp);
end

if strcmpi(prepro_method,'scale') && (NoofSamples >= NoofVariables)
     stdx=std(X.raw);     
     Xecvamodel=X.raw./(ones(size(X.raw,1),1)*stdx);
     % ECVAm.Detail.stdx=stdx;
     ECVAm.Detail.ECVAmodel=sub_ecva(Xecvamodel,X.ClassmarkerNumeric,min_no_of_comp,max_no_of_comp);
else
     ECVAm.Detail.ECVAmodel=sub_ecva(X.raw,X.ClassmarkerNumeric,min_no_of_comp,max_no_of_comp);
end
% ECVAm.Detail.mx=mean(X.raw);

ECVAm.CanonicalWeights=ECVAm.Detail.ECVAmodel.CVA_weights{ECVAm.OptComp};
ECVAm.CanonicalVariates=ECVAm.Detail.ECVAmodel.ProjScores{ECVAm.OptComp};

% Outlier check - residuals
for i=min_no_of_comp:max_no_of_comp
    Wtemp=ECVAm.Detail.ECVAmodel.CVA_weights{i};
    Etransposed=(sub_mc(X.raw) - sub_mc(X.raw)*Wtemp*pinv(Wtemp'*Wtemp)*Wtemp')';
    Q2(:,i)=sum(Etransposed.^2);
end
ECVAm.OutlierDiag.SqRes=Q2;
ECVAm.OutlierDiag.SqResLim=mean(Q2)+3*std(Q2);

% Outlier check - Leverage/Hotelling
for i=min_no_of_comp:max_no_of_comp
    % CVAtemp=ECVAm.Detail.ECVAmodel.ProjScores{i};
    % T2=[];
end
% ECVAm.OutlierDiag.T2=T2;
% ECVAm.OutlierDiag.T2Lim=mean(T2)+3*std(T2);

% Calculate auto-prediction (=fit=cal) and if no validation delete MisClassVal
switch cv_method
    case {'full', 'syst111', 'syst123', 'random', 'manual'}
        ECVAm.MisClassCal = NaN*zeros(1,max_no_of_comp);
        
        Xnoval=X.raw;
        if strcmpi(prepro_method,'scale') && (NoofSamples >= NoofVariables)
           stdx=std(Xnoval);     
           Xnoval=Xnoval./(ones(size(X.raw,1),1)*stdx);
        end
        ECVAall=sub_ecva(Xnoval,X.ClassmarkerNumeric,min_no_of_comp,max_no_of_comp);
        % Mean center calibration objects before projection onto W
        [Xnoval,mx]=sub_mc(Xnoval);
        for j=min_no_of_comp:max_no_of_comp
            PredScores{j}=Xnoval*ECVAall.CVA_weights{j}; % Autoprediction, i.e. fit
            ClassPred{j}=sub_lda(PredScores{j},ECVAall.ProjScores{j},X.ClassmarkerNumeric,prior);
            ECVAm.MisClassCal(j)=length(X.ClassmarkerNumeric)-size(find((ClassPred{j} - X.ClassmarkerNumeric)==0),1);
        end
    case {'none'} % No validation
        ECVAm.MisClassCal=ECVAm.MisClassVal;
        ECVAm.MisClassVal=[];
end

ECVAm.Table=sub_classtable(ECVAm);

% Data compression - scaling back
if NoofSamples < NoofVariables
    ECVAm.CanonicalWeights = V*ECVAm.CanonicalWeights;
    for j=min_no_of_comp:max_no_of_comp
        ECVAm.Detail.ECVAmodel.CVA_weights{j} = V*ECVAm.Detail.ECVAmodel.CVA_weights{j};
    end
end
% Remember to scale back Sw, Sw, St

% Plot
if nargin<=7 || (nargin==8 && plots~=0)
    if X.noofgroups>2 && size(ECVAm.Detail.ECVAmodel.Sb,1)>1
        ecva('cv2d',ECVAm,1,2)
    else
        ecva('cv1d',ECVAm,1)
    end
end

varargout{1}=ECVAm;

function varargout = b_ecvapred(varargin)
%
%  Prediction(s) based on a calculated ECVA model
%
%  Input:
%    Xpred:    a double array containing the samples to be predicted
%    ECVAm:    the output from ecva based on ('model')
%    no_of_comp: [optional] the classification is performed for this number of latent
%                variables in the inner PLS relation. If left out the optimal number
%                as determined by the ECVAcv file is used
%
%  Output:
%    ECVAp: struct array with three fields:
%             .ClassPred: the class predictions
%             .CompUsedForPrediction: the number of components used for prediction
%             .Prior: the prior used in ecva based on ('model')
%
%  ECVAp = ecva('pred',Xpred,ECVAm,no_of_comp)
%    or
%  ECVAp = ecva('pred',Xpred,ECVAm)

Xpred=varargin{1};ECVAm=varargin{2};
switch nargin
    case 3
        no_of_comp=varargin{3};
end

if nargin==2
    no_of_comp=ECVAm.OptComp;
end

if strcmpi(ECVAm.PreproMethod,'scale')
    % The order is important: Should correspond to the order ECVAm.Detail.mx and ECVAm.Detail.stdx are found    
    Xpred=Xpred-(ones(size(Xpred,1),1)*ECVAm.Detail.mx);
    Xpred=Xpred./(ones(size(Xpred,1),1)*ECVAm.Detail.stdx);
else
    Xpred=Xpred-(ones(size(Xpred,1),1)*ECVAm.Detail.mx);
end

PredScores = Xpred*ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp};
PredResiduals = Xpred-Xpred*ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}*pinv(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}'*ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp})*ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}';
PredSqRes = sum(PredResiduals'.^2)';
Bad = PredSqRes>ECVAm.OutlierDiag.SqResLim(no_of_comp);

PredClassNumeric = sub_lda(PredScores,ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp},ECVAm.ClassmarkerNumeric,ECVAm.Prior);

ECVAp.PredClassNumeric=PredClassNumeric;

ECVAp.PredClassOriginal = zeros(length(ECVAp.PredClassNumeric),1);
ClassmarkerOriginalUnique = unique(ECVAm.ClassmarkerOriginal);
ClassmarkerNumericUnique = unique(ECVAm.ClassmarkerNumeric);
for i = 1:size(Xpred,1)
    ixtemp = ECVAp.PredClassNumeric(i) == ClassmarkerNumericUnique;
    ECVAp.PredClassOriginal(i) = ClassmarkerOriginalUnique(ixtemp);
end

ECVAp.CompUsedForPrediction=no_of_comp;
ECVAp.Prior=ECVAm.Prior;
ECVAp.Outliers.PredSqRes=PredSqRes;
ECVAp.Outliers.Bad=Bad;

varargout{1} = ECVAp;

function c1_cv1d(varargin)
%
%  Plots the canonical variates
%
%  Input:
%    ECVAm: results from ecva based on 'model' input
%    cv_a: canonical variate to plot
%    SampleNames [optional]: a char array with variable labels (will be plotted with colors)
%                            if [] is input: only colors will be added to the plot
%                            if 'numbers' is input: colors plus numbers will be added
%                            Only implemented for 2D plots
%    no_of_comp [optional]: number of components to use in inner PLS model
%                         if left out the optimal number estimated by ECVAcv is used
%
%  ecva('cv1d',ECVAm,cv_a,SampleNames,no_of_comp)
%  Examples:
%  ecva('cv1d',ECVAm,1)
%    or
%  ecva('cv1d',ECVAm,1,[],3)
%    or
%  ecva('cv1d',ECVAm,1,SampleNames)

ECVAm=varargin{1};ecv_a=varargin{2};
switch nargin
    case 3
        % SampleNames=varargin{3}; % Not in use for 1d plots(?)
    case 4
        % SampleNames=varargin{3};
        no_of_comp=varargin{4};
end

if nargin<=3
    no_of_comp=ECVAm.OptComp;
end

if ecv_a > length(ECVAm.Prior)-1
    disp(' ')
    stext = sprintf('   Number of canonical variates should be less than the number of groups (%g in this case)',length(ECVAm.Prior));
    disp(stext)
    disp(' ')
    return
end

if length(ecv_a)>1
    disp(' ')
    disp('   This function can only plot one canonical variate at a time')
    disp(' ')
    return
end

Color=['b' 'g' 'r' 'c' 'm' 'y' 'k' 'b' 'g' 'r' 'c' 'm' 'y' 'k']; % Same as MATLAB order
set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
% bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];
set(gcf,'Position',figpos)
bar(1:length(ECVAm.ClassmarkerNumeric),ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp}(:,ecv_a),'FaceColor','w','EdgeColor','w')
hold on
for i = 1:length(unique(ECVAm.ClassmarkerNumeric))
    ixGrp = find(ECVAm.ClassmarkerNumeric==i);
    h(i) = bar(ixGrp,ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp}(ixGrp,ecv_a),Color(i));
    ix(i) = ixGrp(i);
end
hold off
set(gca,'xlim',[0 size(ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp},1)+1])
title('Extended Canonical Variates')

if isnumeric(ECVAm.ClassmarkerOriginal)
    legend(h,int2str(ECVAm.ClassmarkerOriginal(ix')),'Location','NorthEastOutside')
else
    legend(h,ECVAm.ClassmarkerOriginal(ix'),'Location','NorthEastOutside')
end
xlabel('Sample number')
s = sprintf('ECV# %g',ecv_a);
ylabel(s)

function c2_cv2d(varargin)
%
%  Plots the canonical variates
%
%  Input:
%    ECVAm: results from ecva based on 'model' input
%    cv_a: canonical variate to plot on the abscissa
%    cv_b: canonical variate to plot on the ordinate
%    SampleNames [optional]: a char array with variable labels (will be plotted with colors)
%                            if [] is input: only colors will be added to the plot
%                            if 'numbers' is input: colors plus numbers will be added
%                            Only implemented for 2D plots
%    no_of_comp [optional]: number of components to use in inner PLS model
%                         if left out the optimal number estimated by ECVAcv is used
%    NOTE: If two groups are discriminated only one CV exist. Then use 'cv1d'
%
%  ecva('cv2d',ECVAm,cv_a,cv_b,SampleNames,no_of_comp)
%  Examples:
%  ecva('cv2d',ECVAm,1,2)
%    or
%  ecva('cv2d',ECVAm,1,2,[],3)
%    or
%  ecva('cv2d',ECVAm,1,2,SampleNames)

ECVAm=varargin{1};cv_a=varargin{2};cv_b=varargin{3};
switch nargin
    case 4
        SampleNames=varargin{4};
    case 5
        SampleNames=varargin{4};
        no_of_comp=varargin{5};
end

if nargin<=4
    no_of_comp=ECVAm.OptComp;
end

if length(ECVAm.Prior)==2 % Two groups
    disp(' ')
    disp('   Use ecva(''cv1d'') to show plot for two group models')
    disp(' ')
    return
end

if cv_a > size(ECVAm.CanonicalWeights,2) || cv_b > size(ECVAm.CanonicalWeights,2)
  disp(' ')
  disp('   The highest canonical variate number should be lower than the number of groups')
  disp(' ')
  return
end

if nargin <=3 || isempty(SampleNames)
    sub_scatter(ECVAm,cv_a,cv_b,[],no_of_comp);
else
    sub_scatter(ECVAm,cv_a,cv_b,SampleNames,no_of_comp);
end

function d_cw1d(varargin)
%
%  Line plot of the weights
%
%  Input:
%    ECVAm:  results from ecva with 'model' as input
%    weight: insert either a number, e.g. 1, or several numbers, e.g. 1:3 or [1:2 5]
%    no_of_comp: [optional] number of components to use in inner PLS model
%                if left out the optimal number estimated by ECVAcv is used
%
%  ecva('cw1d',ECVAm,weight,no_of_comp)
%  Examples:
%  ecva('cw1d',ECVAm,1)
%    or
%  ecva('cw1d',ECVAm,2:3)

ECVAm=varargin{1};weight=varargin{2};
switch nargin
    case 3
        no_of_comp=varargin{3};
end

if max(weight) > size(ECVAm.CanonicalWeights,2)
    disp(' ')
    disp('   The highest canonical weight number should be lower than the number of groups')
    disp(' ')
    return
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
% bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];
set(gcf,'Position',figpos)

if nargin==2
    no_of_comp=ECVAm.OptComp;
end

if size(ECVAm.CanonicalWeights,2)>1 && no_of_comp<size(ECVAm.CanonicalWeights,2)
  % CHECK ABOVE: Seems to work with 4 comp for fluo data with 6 groups??
  no_of_comp=size(ECVAm.CanonicalWeights,2);
  s=sprintf('No_of_lv is set to ''number of groups minus one: %g',no_of_comp); disp(s);
end

plot(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,weight))
% bar(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,weight))
title('Extended Canonical Weights')
axis('tight')
xlabel('Variable number')

% Preallocating
strcell=cell(size(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,weight),2),1);

for i=1:size(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,weight),2)
    str=['ECW# ' int2str(weight(i))];
    strcell{i}=str;
end

legend(strcell);
sub_addcross;

function e_cw2d(varargin)
%
%  2D plot of the weights
%
%  Input:
%    ECVAm: results from ECVA
%    w_a: the first weight vector
%    w_b: the second weight vector
%    VarNames: [optional] a char array with variable labels, use []
%              if empty then numbers will be added
%    no_of_comp: [optional] number of components to use in inner PLS model
%              if left out the optimal number estimated by ECVAcv is used
%    NOTE: this plot is not relevant for two group data sets
%
%  ecva('cw2d',w_a,w_b,VarNames,no_of_comp)
%  Examples:
%  ecva('cw2d',ECVAm,1,2,[],5)
%    or
%  ecva('cw2d',ECVAm,1,2)

ECVAm=varargin{1};w_a=varargin{2};w_b=varargin{3};
switch nargin
    case 4
        VarNames=varargin{4};
    case 5
        VarNames=varargin{4};
        no_of_comp=varargin{5};
end

if size(ECVAm.CanonicalWeights,2)==1
  disp(' ')
  disp('   Plot not relevant for two group data')
  disp(' ')
  return
end

if w_a > size(ECVAm.CanonicalWeights,2) || w_b > size(ECVAm.CanonicalWeights,2)
  disp(' ')
  disp('   The highest canonical weight number should be lower than the number of groups')
  disp(' ')
  return
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
% bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];
set(gcf,'Position',figpos)

if nargin<5
  no_of_comp=ECVAm.OptComp;
end

if nargin==3
   VarNames=[];
end

if size(ECVAm.CanonicalWeights,2)>1 && no_of_comp<size(ECVAm.CanonicalWeights,2)
  % NOTE: Check above; seems to work with 4 comp for fluo data with 6 groups?
  no_of_comp=size(ECVAm.CanonicalWeights,2);
  s=sprintf('No_of_lv is set to ''number of groups minus one'' = %g',no_of_comp); disp(s);
end

n=length(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_a));
plot(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_a),ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_b),'w')
if isempty(VarNames)
    text(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_a),ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_b),int2str((1:n)'))
else
    text(ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_a),ECVAm.Detail.ECVAmodel.CVA_weights{no_of_comp}(:,w_b),VarNames)    
end
title('Extended Canonical Weights')
s=sprintf('ECW# %g',w_a);
xlabel(s)
s=sprintf('ECW# %g',w_b);
ylabel(s)
sub_addcross;

function varargout=f_iecva(varargin)

%  Calculates interval ECVA models (a la iPLS)
%
%  Input:
%    X: the independent variables
%    classmarker: the group variable,
%    no_of_comp: the maximum number of PLS components in the inner relation
%    prepro_method (for X): 'none' (e.g. spectra) and 'scale' (1/Sdev, e.g. chemical data)
%    intervals: the number of intervals
%      if intervals is a row vector divisions are made based on the elements
%      [startint1 endint1 startint2 endint2 startint3 endint3], see an example in manint
%    xaxislabels: for labelling of e.g. spectroscopic data, if not available type []
%    val_method is 'full', 'syst111', 'syst123', 'random', 'manual' or 'none';
%    segments (segments = total number of samples corresponds to full cv)
%      if intervals is a cell array cross validation is performed according to this array
%    prior: as in ECVA.
%    plots: if set to 0 plots and text are suppressed,
%           otherwise set to a number different from 0 (e.g. 1)
%    minusglobal: if set to 0 does not calculate global model
%                 (could be necessary if memory problems occur).
%
%  Output:
%    iECVAm: a structured array containing all model information
%
%  iECVAm=ecva('iecva',X,classmarker,no_of_comp,prepro_method,intervals,xaxislabels,val_method,segments,prior,plots,minusglobal);
%
%  Example:
%  iECVAm=ecva('iecva',X,classmarker,15,'none',20,xaxis,'syst123',5);
%

X=varargin{1};classmarker=varargin{2};no_of_comp=varargin{3};prepro_method=varargin{4};intervals=varargin{5};
xaxislabels=varargin{6};val_method=varargin{7};
prior=[];
plots=1; % Plots as output [default]
minusglobal=1; % Calculates global model [default]
switch nargin
    case 8
        segments=varargin{8};
    case 9
        segments=varargin{8};
        prior=varargin{9};
    case 10
        segments=varargin{8};
        prior=varargin{9};
        plots=varargin{10};
    case 11
        segments=varargin{8};
        prior=varargin{9};
        plots=varargin{10};
        minusglobal=varargin{11};
end

% Error checks
if ~ismember(val_method,{'test', 'full', 'syst123', 'syst111', 'random', 'manual'})
    disp(' Not allowed validation method')
    return
end

if ~ismember(prepro_method,{'none', 'scale'})
    disp(' Not allowed preprocessing method')
    return
end

if strcmpi(val_method,'manual') && ~iscell(segments) %strcmpi ignores case
    disp('You need to specify the manual segments in a cell array, see script_manseg')
    return
end

if strcmpi(val_method,'manual') && iscell(segments)
   Nsamples=0;
   for j=1:max(size(segments))
      Nsamples=Nsamples + max(size(segments{j}));
   end
   if size(X,1)~=Nsamples
      disp('The number of samples in X does not correspond to the total number of samples in segments')
      return
   end
end
% End error checks

iECVAm.type='iECVA';
iECVAm.rawX=X;
iECVAm.classmarker=classmarker;
iECVAm.no_of_comp=no_of_comp;
iECVAm.prepro_method=prepro_method;
iECVAm.xaxislabels=xaxislabels; % Final use of xaxislabels in this file; i.e. if reversed axes are present this should be taken care of in iplot
iECVAm.val_method=val_method;
iECVAm.prior=prior;
iECVAm.minusglobal=minusglobal;

[n,m]=size(X);
if strcmpi(iECVAm.val_method,'full') || nargin <8
   iECVAm.segments=n;
else
   iECVAm.segments=segments;
end

mint=size(intervals,2);
if mint > 1
    iECVAm.allint=[(1:round(mint/2)+1)' [intervals(1:2:mint)';1] [intervals(2:2:mint)';m]];
    iECVAm.intervals=round(mint/2);
    iECVAm.intervalsequi=0;
else
    iECVAm.intervals=intervals;
    vars_left_over=mod(m,intervals);
    N=fix(m/intervals);
    % Distributes vars_left_over in the first "vars_left_over" intervals
    startint=[(1:(N+1):(vars_left_over-1)*(N+1)+1)'; ((vars_left_over-1)*(N+1)+1+1+N:N:m)'];
    endint=[startint(2:intervals)-1; m];
    iECVAm.allint=[(1:intervals+1)' [startint;1] [endint;m]];
    iECVAm.intervalsequi=1;
end

% Local calibration
for i=1:size(iECVAm.allint,1)
   if i<size(iECVAm.allint,1)
       if plots~=0
          home, s = sprintf('Working on interval no. %g of %g...',i,size(iECVAm.allint,1)-1); disp(s)
       end
      % Input 0 supresses plots
      iECVAm.ECVAmodel{i}=ecva('model',X(:,iECVAm.allint(i,2):iECVAm.allint(i,3)),classmarker,iECVAm.no_of_comp,iECVAm.prepro_method,iECVAm.val_method,iECVAm.segments,prior,0);
   else
      if minusglobal==0
         if plots~=0
             home, disp('Uses last interval model as global model for technical reasons ...')
         end
         iECVAm.ECVAmodel{i}=iECVAm.ECVAmodel{i-1};
      else
         if plots~=0
             home, disp('Calculating global model ...')
         end
         iECVAm.ECVAmodel{i}=ecva('model',X(:,iECVAm.allint(i,2):iECVAm.allint(i,3)),classmarker,iECVAm.no_of_comp,iECVAm.prepro_method,iECVAm.val_method,iECVAm.segments,prior,0);
      end
   end
end

varargout{1}=iECVAm;

% Plot results
if plots~=0
    if iECVAm.intervalsequi==0
        ecva('iplot',iECVAm,'varlabel');
    else
        ecva('iplot',iECVAm);
    end
end

function g_iplot(varargin)

%  Bar plot of results from interval ECVA analysis
%
%  Input:
%    iECVAm: the output from ecva('iecva')
%    labeltype [optional]: designates whether you want:
%               interval number ('intlabel'), [default value]
%               variable number ('varlabel')
%               wavelength number ('wavlabel')
%    optimal_lv_global [optional]: the number of PLS components chosen for full spectrum model;
%               if not given or given by [] the minimum classification error is chosen
%    max_yaxis [optional]: control ordinate scaling in the iECVA plot
%    plottype: [optional]
%               'Cum' (default)
%               'n' (where is a locked number of inner PLS components)
%
%  ecva('iplot',iECVAm,labeltype,optimal_lv_global,max_yaxis,plottype);
%    or
%  ecva('iplot',iECVAm,labeltype);
%  Examples:
%  ecva('iplot',iECVAm,'intlabel',[],[],'Cum');
%    or
%  ecva('iplot',iECVAm,'intlabel');
%

iECVAm=varargin{1};labeltype='intlabel';
switch nargin
    case 2
        labeltype=varargin{2};
    case 3
        labeltype=varargin{2};
        optimal_lv_global=varargin{3};
    case 4
        labeltype=varargin{2};
        optimal_lv_global=varargin{3};
        max_yaxis=varargin{4};
    case 5
        labeltype=varargin{2};
        optimal_lv_global=varargin{3};
        max_yaxis=varargin{4};
        plottype=varargin{5};        
end

% Error checks
if ~isstruct(iECVAm)
    disp('Model input should be a structure array (output from iecva.m)')
    return
end

if ~strcmpi(iECVAm.type,'iECVA')
    disp('The model input is not from an iECVA model')
    return
end

plottype_cell=cell(1+iECVAm.no_of_comp,1);
for i=1:iECVAm.no_of_comp
    plottype_cell{i+1}=num2str(i);
end
plottype_cell{1}='Cum';

if nargin<5
    plottype='Cum';    
end

if nargin==5 && ~ismember(plottype,plottype_cell)
    disp('Not legal plottype, use ''Cum'' or different number of PLS components')
    return
end

if ~ismember(labeltype,{'intlabel','varlabel','wavlabel'})
    disp('Not legal labeltype, use ''intlabel'',''varlabel'', or ''wavlabel''')
    return
end

if iECVAm.intervalsequi==0 && strcmpi(labeltype,'intlabel')
   disp(' ')
   disp(' Manually chosen intervals are not correctly plotted with ''intlabel''')
   disp(' so please use ''varlabel'' or ''wavlabel'' as labeltype')
   disp(' ')
end
% End error checks

Xmean = mean(iECVAm.rawX); % Mean spectrum
if min(Xmean)<0
   Xmean=Xmean+(-min(Xmean)); % To make all intensities positive
end
m = size(iECVAm.rawX,2);

No_Int = iECVAm.intervals;
if strcmpi(labeltype,'intlabel')
   Xtext = sprintf('Interval number');
   Xint =[iECVAm.allint(1:No_Int,1)-0.5 iECVAm.allint(1:No_Int,1)-0.5 iECVAm.allint(1:No_Int,1)+0.5 iECVAm.allint(1:No_Int,1)+0.5]';
   NumberofTicksandWhere=mean(Xint(2:3,:));
elseif strcmpi(labeltype,'wavlabel')
   if isempty(iECVAm.xaxislabels)
      disp('You must define wavelength/wavenumber labels')
      return
   end
   Xtext = sprintf('Wavelength/Wavenumber');
   a = iECVAm.allint(1:No_Int,2);
   b = iECVAm.allint(1:No_Int,3);

   % To reverse wavenumber axis before plotting; will be reversed back when the
   % final plot is made
   NewAxisLabels=iECVAm.xaxislabels; % Important; original axislabels are used in the last three lines of the program
   if NewAxisLabels(1)>NewAxisLabels(2)
       if size(NewAxisLabels,1)==1
           NewAxisLabels=fliplr(NewAxisLabels);
       elseif size(NewAxisLabels,2)==1
           NewAxisLabels=flipud(NewAxisLabels);
       end
   end
   
   Xint = [NewAxisLabels(a)' NewAxisLabels(a)' NewAxisLabels(b)' NewAxisLabels(b)']';
   NumberofTicksandWhere=[Xint(2,:) Xint(3,end)];

elseif strcmpi(labeltype,'varlabel')
   Xtext = sprintf('Variable number');
   Xint  = [iECVAm.allint(1:No_Int,2) iECVAm.allint(1:No_Int,2) iECVAm.allint(1:No_Int,3) iECVAm.allint(1:No_Int,3)]';
   NumberofTicksandWhere=[Xint(2,:) Xint(3,end)];
end

switch plottype
    case {'Cum'}
		for i=1:iECVAm.intervals
            min_ix(i)=iECVAm.ECVAmodel{i}.OptComp;
            minMisClass(i)=iECVAm.ECVAmodel{i}.MisClassVal(min_ix(i));
        end
    case plottype_cell(2:end)
        for i=1:iECVAm.intervals
            min_ix(i)=str2num(plottype); % Same number of components for all intervals
            minMisClass(i)=iECVAm.ECVAmodel{i}.MisClassVal(min_ix(i));
        end
    otherwise
        disp('Not valid plot type')
end % switch

if nargin<=2
  optimal_lv_global=iECVAm.ECVAmodel{end}.OptComp;
end

if nargin >2 && isempty(optimal_lv_global)
  optimal_lv_global=iECVAm.ECVAmodel{end}.OptComp;
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
% bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];
figure(1)
set(1,'Position',figpos)
Response = [zeros(No_Int,1) minMisClass(1:No_Int)' minMisClass(1:No_Int)' zeros(No_Int,1)]';

% Cumulated plots
if nargin <5 || ismember(plottype,plottype_cell([1 2:end]))
    if strcmpi(labeltype,'wavlabel') && (iECVAm.xaxislabels(1)>iECVAm.xaxislabels(2))
        fill(flipud(Xint(:)),Response(:),[0.75 0.75 0.75])
    else
        fill(Xint(:),Response(:),[0.75 0.75 0.75]) % Note: substitute [0.75 0.75 0.75] for 'c'
    end
else
    % Nothing
end

if strcmpi(iECVAm.val_method,'test')
    if iECVAm.minusglobal~=0
        plottitle = sprintf('Dotted line is no of misclassifications (%g for %g LV''s) for global model / Italic numbers are optimal LVs in interval model',iECVAm.ECVAmodel{end}.MisClassVal(optimal_lv_global),optimal_lv_global);
    else
        plottitle = sprintf('No global model calculated / Italic numbers are optimal LVs in interval model');
    end
    ylabel('No of misclassifications / test set','FontSize',10)    
elseif strcmpi(iECVAm.val_method,'none')
    if iECVAm.minusglobal~=0
        plottitle = sprintf('Dotted line is no of misclassifications (%g for %g LV''s) for global model / Italic numbers are optimal LVs in interval model',iECVAm.ECVAmodel{end}.MisClassVal(optimal_lv_global),optimal_lv_global);
    else
        plottitle = sprintf('No global model calculated / Italic numbers are optimal LVs in interval model');
    end
    ylabel('No of misclassifications / no validation','FontSize',10)
else
    if iECVAm.minusglobal~=0
        plottitle = sprintf('Dotted line is no of misclassifications (%g for %g LV''s) for global model / Italic numbers are optimal LVs in interval model',iECVAm.ECVAmodel{end}.MisClassVal(optimal_lv_global),optimal_lv_global);
    else
        plottitle = sprintf('No global model calculated / Italic numbers are optimal LVs in interval model');
    end
    ylabel('No of misclassifications / cross validation','FontSize',10)
end
title(plottitle,'FontSize',10,'FontWeight','Bold')
xlabel(Xtext)

hold on
    axis tight;
    if iECVAm.minusglobal~=0
        sub_horzline(iECVAm.ECVAmodel{end}.MisClassVal(optimal_lv_global),':k')
    end
    actualaxis=axis;
    if nargin >= 4 && ~isempty(max_yaxis)
        axis([actualaxis(1) actualaxis(2) actualaxis(3) max_yaxis]);
        actualaxis(4)=max_yaxis;
    end
    Xaxis = linspace(actualaxis(1),actualaxis(2),m);
    if strcmpi(labeltype,'wavlabel') && (iECVAm.xaxislabels(1)>iECVAm.xaxislabels(2))
        plot(fliplr(Xaxis),Xmean./max(Xmean)*actualaxis(4),'-k') % Scaled spectrum
    else
        plot(Xaxis,Xmean./max(Xmean)*actualaxis(4),'-k') % Scaled spectrum
    end
    set(gca,'XTick',NumberofTicksandWhere)
    for i=1:iECVAm.intervals
      	if strcmpi(labeltype,'wavlabel') && (iECVAm.xaxislabels(1)>iECVAm.xaxislabels(2))
            text(mean(Xint(2:3,i)),0.03*(actualaxis(4)-actualaxis(3))+actualaxis(3),int2str(min_ix(iECVAm.intervals-(i-1))),'Color','k','FontAngle','italic');
        else
            text(mean(Xint(2:3,i)),0.03*(actualaxis(4)-actualaxis(3))+actualaxis(3),int2str(min_ix(i)),'Color','k','FontAngle','italic');
        end
    end
hold off

% To reverse in case of reversed wavenumber axis
if strcmpi(labeltype,'wavlabel') && (iECVAm.xaxislabels(1)>iECVAm.xaxislabels(2))
    set(gca,'XDir','reverse');
end

function h_icanonvar(varargin)
%
%  2D (or 1D if two groups) plots of ECVa versus ECVb for all intervals
%
%  Input:
%    iECVAm: the output from ecva('iecva')
%    ecv_a and ecv_b designates the combination of canonical variates to use in the plot
%           If only two groups: just use the same values for ecv_a and ecv_b
%    SampleNames [optional]: char, use [] if not available
%    no_of_plot_pr_figure [optional]: must be 2, 4, 6 or 8 (optional, default is 6)
%
%  ecva('icanonvar',iECVAm,ecv_a,ecv_b,samplenames,no_of_plots_pr_figure);
%  Example:
%  ecva('icanonvar',iECVAm,1,2,[]);
%

iECVAm=varargin{1};ecv_a=varargin{2};ecv_b=varargin{3};
SampleNames=[];
no_of_plots_pr_figure=6;
switch nargin
    case 4
        SampleNames=varargin{4};
    case 5
        SampleNames=varargin{4};
        no_of_plots_pr_figure=varargin{5};
end

% Plot score plot for global model
figure(1)
if length(unique(iECVAm.classmarker))==2
    if ecv_a~=1
        ecv_a=1;
        disp(' ')
        disp('   ecv_a is set to 1 (only two groups are present)')
        disp(' ')        
    end
    ecva('cv1d',iECVAm.ECVAmodel{end},ecv_a,SampleNames)
else
    ecva('cv2d',iECVAm.ECVAmodel{end},ecv_a,ecv_b,SampleNames)
end

if isempty(iECVAm.xaxislabels)==1
   title(sprintf('Global model, Variables %g-%g',iECVAm.allint(1,2),iECVAm.allint(size(iECVAm.allint,1),3)));
else
   wavstart=iECVAm.xaxislabels(iECVAm.allint(1,2));
   wavend=iECVAm.xaxislabels(iECVAm.allint(size(iECVAm.allint,1),3));
   title(sprintf('Global model, Var. %g-%g, Wav. %g-%g',iECVAm.allint(1,2),iECVAm.allint(size(iECVAm.allint,1),3),wavstart,wavend));
end

% Plot score plots for local models
No_of_figures=ceil(iECVAm.intervals/no_of_plots_pr_figure);
count=0;
for i=1:No_of_figures
   figure(i+1)
   for j=1:no_of_plots_pr_figure
      count=count+1;
      if count>iECVAm.intervals, break, end
      if no_of_plots_pr_figure==1
          % One plot in one Figure
      elseif no_of_plots_pr_figure==2
          subplot(1,2,j)
      elseif no_of_plots_pr_figure==4
          subplot(2,2,j)
      elseif no_of_plots_pr_figure==8
          subplot(4,2,j)          
      else   %no_of_plots_pr_figure==6
          subplot(3,2,j)
      end
      if length(unique(iECVAm.classmarker))==2
          ecva('cv1d',iECVAm.ECVAmodel{count},ecv_a,SampleNames)
      else
          ecva('cv2d',iECVAm.ECVAmodel{count},ecv_a,ecv_b,SampleNames)
      end
      if isempty(iECVAm.xaxislabels)==1
        title(sprintf('Int%g, Variables %g-%g',count,iECVAm.allint(count,2),iECVAm.allint(count,3)));
      else
        wavstart=iECVAm.xaxislabels(iECVAm.allint(:,2));
        wavend=iECVAm.xaxislabels(iECVAm.allint(:,3));
        title(sprintf('Int%g, Var. %g-%g, Wav. %g-%g',count,iECVAm.allint(count,2),iECVAm.allint(count,3),wavstart(count),wavend(count)));
      end
   end
end

function i_intervals(varargin)

%
%  Prints (on screen) the intervals with variable number and/or wavelength label
%  Input is an interval model (the output from ecva('iecva'))
%
%  ecva('intervals',iECVAm);

iECVAm=varargin{1};

disp(' ')
if isempty(iECVAm.xaxislabels)
  disp('     Int.no     Start        End      No. vars.')
  number_of_vars=iECVAm.allint(1:end-1,3)-iECVAm.allint(1:end-1,2)+1;
  table=[iECVAm.allint(1:end-1,:) number_of_vars];
  fprintf(1,'%10.0f %10.0f %10.0f %10.0f\n',table')
else
  disp('    Interval  Start var.   End var.  Start wav.  End wav.  Number of vars.')
  number_of_vars=iECVAm.allint(1:end-1,3)-iECVAm.allint(1:end-1,2)+1;
  table=[iECVAm.allint(1:end-1,:) iECVAm.xaxislabels(iECVAm.allint(1:end-1,2))' iECVAm.xaxislabels(iECVAm.allint(1:end-1,3))' number_of_vars];
  fprintf(1,'%10.0f %10.0f %10.0f %10.0f %10.0f %10.0f\n',table')
end

function varargout = k_randtest(varargin)
%
%  Calculates an ECVA model on random X data of the same dimension as the
%    input and the same settings as used in the ECVA model of the real data
%
%  Input:
%    X: a double array containing the independent variables
%    classmarker (column vector): contains the group information (you can use either numbers or strings).
%    no_of_comp: the number of PLS components in the inner CVA relation
%              NOTE: choose a large number since this PLS model includes a variables x variables matrix
%              The optimal number will be found by cross validation (if this option is selected)
%    prepro_method: 'none' (e.g. NIR spectra) or 'scale' (1/SDev, e.g. chemical data)
%    cv_method: 'full' (leave-one-out), 'syst111' (cont. blocks), 'syst123' (Venetian blinds), 'random', 'manual' or 'none'
%               if 'none' is selected, insert [] in segments (or only use five inputs)
%    segments: segments in cross validation
%              Using the number of samples in X corresponds to full cv ('full' overrules segments)
%              if cv_methods is 'manual' use script_manseg to form the segments,
%              then insert the cell array formed (called manualseg) instead of the number of segments
%    prior: [optional], if set to [] equal group probability will be used
%    plots: [optional], if set to 0 no plots & text will be output, otherwise set to 1
%
%  Output:
%    ECVAm: a structured array containing all model and validation information
%             ECVAm.MisClassVal contains the number of misclassifications
%               (validated) as a function of the number of PLS components
%             ECVAm.ClassPredOptComp to see the class predictions for the
%               optimal number of PLS components
%
%  ECVAm = ecva('randtest',X,classmarker,no_of_comp,prepro_method,cv_method,segments,prior,plots);
%  Examples: 
%  ECVAm = ecva('randtest',X,classmarker,10,'none','syst123',15,[],0);
%    or
%  ECVAm = ecva('randtest',X,classmarker,10,'none','syst123',15,[0.25 0.25 0.50]);
%    or
%  ECVAm = ecva('randtest',X,classmarker,10,'scale','syst123',15);

varargin{1} = randn(size(varargin{1}));
ECVAm = ecva('model',varargin{:});
ECVAm.RandTest.CalErrorAverageOverPLSComp = mean( ECVAm.MisClassCal(size(unique(varargin{2}),1):end)) / size(varargin{1},1);
ECVAm.RandTest.ValErrorAverageOverPLSComp = mean( ECVAm.MisClassVal(size(unique(varargin{2}),1):end)) / size(varargin{1},1);
ECVAm.RandTest.ExpectedErrorIfRandom = 1 - (1/size(unique(varargin{2}),1)); % No prior or number of samples taken into account!
varargout{:} = ECVAm;

function ECVAsub = sub_ecva(X,classmarker,min_no_of_comp,max_no_of_comp)

%  sub_ecva, used by ECVAcv
%  X: double array
%  classmarker: vector with class information
%  min_no_of_comp: minimum number of PLS components to extract in inner model
%  max_no_of_comp: maximum number of PLS components to extract in inner model
%  NOTE: if you want to work on scaled (1/Sdev) data use ECVAcv with the parameter 'scale'
%        or scale the data before application of this file (DO NOT autoscale, only scale)!

X = sub_makestruct(X,classmarker);

% Initialize and calculate Swithin
ECVAsub.Sw = zeros(size(X.raw,2),size(X.raw,2));
for i = 1:X.noofgroups
    ECVAsub.Sw = ECVAsub.Sw + cov(X.raw(X.ix{i},:))*(X.N{i}-1);
end

% Calculate Stotal
ECVAsub.St = cov(X.raw)*(X.NoofSamples-1);

% Calculate Sbetween
ECVAsub.Sb = ECVAsub.St - ECVAsub.Sw;

m=mean(X.raw);
for i = 1:X.noofgroups
    Means{i}=mean(X.raw(X.ix{i},:));
    y(:,i)=(Means{i}-m)';
end

% Special case of two groups
if X.noofgroups==2
    m1=mean(X.raw(X.ix{1},:));
    m2=mean(X.raw(X.ix{2},:));
    y=(m1-m2)';
end

% Calculate the weights through PLS models
% Note: with e.g. six groups, six y-variables (and not five) are used to obtain an
%       equal weighting of directions in PLS2
% It seems that NO mean centering gives sligthly better classification

ECVAsub.Bw = sub_simpls(ECVAsub.Sw,y,max_no_of_comp);
CVA_weights = ECVAsub.Bw;

% At this point the number of CVAs will be equal to number of groups except
%    if the input number of groups is 2, then the number of CVAs is 1 and no
%    elimination is necessary

% Order the weights according to optimization criterion and eliminate the poorest
if X.noofgroups >2 && X.noofgroups<=size(X.raw,2)
    for i=min_no_of_comp:max_no_of_comp
        for j=1:X.noofgroups
           w_temp=CVA_weights{i}(:,j);
           OptCrit(j) = det(w_temp'*ECVAsub.Sb*w_temp)/det(w_temp'*ECVAsub.Sw*w_temp);
        end
        [~,sort_index]=sort(OptCrit,'descend');
        CVA_weights{i}=CVA_weights{i}(:,sort_index);
        % Flip sign so that the average element is positive; this is done in manova1 in Statistics Toolbox
        switch_index = (sum(CVA_weights{i}) < 0);
        CVA_weights{i}(:,switch_index) = -CVA_weights{i}(:,switch_index);
        % Leave out the poorest performing
        CVA_weights{i}(:,end)=[];
    end
elseif X.noofgroups>size(X.raw,2) % For number of groups larger than number of X variables
    for j=min_no_of_comp:max_no_of_comp
        CVA_weights{j}=CVA_weights{j}(:,1:j); % Since rank of b's for the PLS2 models is equal to number of PLS components
    end
end

% Collect, rotate and normalize weights ("loadings") for each component
for j=min_no_of_comp:max_no_of_comp
    % Rotate weights to wi*Swithin*wj=0
    [eigvecs,~]=eig(CVA_weights{j}'*ECVAsub.Sw*CVA_weights{j});
    CVA_weights{j}=CVA_weights{j}*eigvecs;
    % Normalize weights to wi*Swithin*wi=1 (have to be after rotation)
    tmpMat = CVA_weights{j}'*ECVAsub.Sw*CVA_weights{j};
    k=(ones(size(CVA_weights{j},2),1))./diag(tmpMat);
    ECVAsub.CVA_weights{j}=CVA_weights{j}.*(ones(size(CVA_weights{j},1),1)*sqrt(k)');
end

% Calculate canonical variates - the "scores"
% Mean centering of X.raw is performed as in the standard CVA (manoval)
for j=min_no_of_comp:max_no_of_comp
    ECVAsub.ProjScores{j}=sub_mc(X.raw)*ECVAsub.CVA_weights{j};
    % ECVAsub.ProjScores{j}=X.raw*ECVAsub.CVA_weights{j};
end

% Calculate optimization criterion
for j=min_no_of_comp:max_no_of_comp
    w_temp=ECVAsub.CVA_weights{j};
   ECVAsub.V(j) = det(w_temp'*ECVAsub.Sb*w_temp)/det(w_temp'*ECVAsub.Sw*w_temp);
end

%trSwithin = trace(Swithin); trSbetween = trace(Sbetween); trStotal = trace(Stotal); Traces = [trSwithin trSbetween trStotal];
%Traces = [trSwithin/trStotal trSbetween/trStotal];

function sub_horzline(ordinate,linetype_color)

V = axis;
plot(V(1:2),[1 1]*ordinate,linetype_color);

function [new_group,posterior] = sub_lda(Xnew,Xtrain,classmarker,prior)

%  SUB_LDA Linear Discriminant Analysis for ECVA function
%
%  Input:
%  Xnew: double array containing the new samples to be classified
%  Xtrain: double array containing the calibration samples
%  classmarker (column vector): contains the group information (you can use all sorts of
%           groupnames: numbers, strings etc.). The code locates the unique
%           values and makes groups accordingly with numbers starting from 1.
%  prior [optional]: prior probability for each group; e.g. [0.25 0.25 0.50]
%                    set to [] or omit if equal probabilities are to be used
%
%  Output:
%  new_group: the groups assigned to the samples in Xnew
%
%  new_group = sub_lda(Xnew,Xtrain,classmarker,prior);
%  Example: 
%  new_group = sub_lda(Xnew,Xtrain,classmarker,[0.25 0.25 0.5]);

n=size(Xtrain,1);

if isnumeric(classmarker)
    classmarkerUnique=unique(classmarker);
	ix_class=NaN*ones(n,1);
	for j=1:length(classmarkerUnique)
        for i=1:n
            ixtemp(i)=eq(classmarker(i,:),classmarkerUnique(j,:));
        end
        ix{j}=find(ixtemp==1)';
        N{j}=length(ix{j});
        ix_class(ix{j})=j;
	end
	classmarker=ix_class;
else
	classmarkerUnique=unique(classmarker,'rows');
	ix_class=NaN*ones(n,1);
	for j=1:length(classmarkerUnique)
        for i=1:n
            ixtemp(i)=strcmpi(classmarker(i,:),classmarkerUnique(j,:));
        end
        ix{j}=find(ixtemp==1)';
        N{j}=length(ix{j});
        ix_class(ix{j})=j;
	end
	classmarker=ix_class;
end

noofgroups = length(classmarkerUnique);

if nargin <= 3 || isempty(prior)
    prior = ones(1,noofgroups-1)/noofgroups;
    prior = [prior 1-sum(prior)];
end

if nargin == 4
    if abs(sum(prior)-1) > 100*eps
    % OLD: if sum(prior) ~= 1
        disp('The priors do not sum to one');
    end
end

% Initialize and calculate Swithin
Swithin = zeros(size(Xtrain,2),size(Xtrain,2));
% for i = 1:noofgroups
%     Stemp = cov(Xtrain(ix{i},:))*(N{i}-1);  
%     Swithin = Swithin + Stemp;
%     x_average(i,:)=mean(Xtrain(ix{i},:));
% end
% NEW CODE - RB
for i = 1:noofgroups
    %Stemp = cov(Xtrain(ix{i},:))*(N{i}-1);  
    %Stemp = cov(Xtrain(ix{i},:));%*(N{i}-1);  
    x_average(i,:)=mean(Xtrain(ix{i},:));
    Xm = Xtrain(ix{i},:)-ones(N{i},1)*x_average(i,:);
    Stemp = Xm'*Xm; 
    Swithin = Swithin + Stemp;
end

Swithin = Swithin/(size(Xtrain,1)-noofgroups); % Added by RB
L=zeros(size(Xnew,1),noofgroups);
for i = 1:size(Xnew,1)
    for j = 1:noofgroups
        L(i,j) = log(prior(j)) - 0.5*(Xnew(i,:)-x_average(j,:))*pinv(Swithin)*(Xnew(i,:)-x_average(j,:))' + log(det(Swithin));
    end
end

% Added by RB
[~,new_group]=max(L,[],2);
maxL = max(L,[],2);
P = exp(L - repmat(maxL,1,size(L,2)));
P = exp(L);
sumP = sum(P,2);
posterior = P ./ repmat(sumP,1,size(L,2));

function Xstruct=sub_makestruct(X,classindex)

% Makes X structured array from raw data and class input
% Xstruct=struct('ClassmarkerOriginal',[],'NoofSamples',[],'OptComp',[],'ClassPredOptComp',[],...
%                    'CanonicalWeights',[],'CanonicalVariates',[],'ValidationMethod',[],'PreproMethod',[]);

if isnumeric(classindex) && size(classindex,1)==1
    classindex = classindex';
end

if iscell(classindex)
    classindex = char(classindex);
end

Xstruct.NoofSamples=size(classindex,1);

Xstruct.ClassmarkerOriginal=classindex;

if isnumeric(classindex)
    Xstruct.classmarker_string = int2str(classindex);
    classindexUnique = unique(classindex);
	ix_class = NaN*ones(Xstruct.NoofSamples,1);
	for j = 1:length(classindexUnique)
        for i = 1:Xstruct.NoofSamples
            ixtemp(i) = eq(classindex(i,:),classindexUnique(j,:));
        end
        ix{j} = find(ixtemp == 1)';
        N{j} = length(ix{j});
        ix_class(ix{j}) = j;
	end
	classindex = ix_class;
    Xstruct.noofgroups=length(classindexUnique);
else
	classindexUnique=unique(classindex,'rows');
    ix_class=NaN*ones(Xstruct.NoofSamples,1);
	for j=1:size(classindexUnique,1)
        for i=1:Xstruct.NoofSamples
            ixtemp(i)=strcmpi(classindex(i,:),classindexUnique(j,:));
        end
        ix{j}=find(ixtemp==1)';
        N{j}=length(ix{j});
        ix_class(ix{j})=j;
	end
	classindex=ix_class;
    Xstruct.noofgroups=size(classindexUnique,1);
end
if ( length(unique(classindex)) ~= length(classindexUnique))
    classindex;
end
Xstruct.N=N;
Xstruct.raw=X;
Xstruct.ix=ix;
Xstruct.ClassmarkerNumeric = classindex;

function [Xmean,meanX] = sub_mc(X)
n = size(X,1);
meanX = mean(X);
Xmean = (X-meanX(ones(n,1),:));

function sub_scatter(ECVAm,cv_a,cv_b,SampleNames,no_of_comp)

if isnumeric(SampleNames)
    SampleNames=int2str(SampleNames);
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
% bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.1*ScrLength 0.15*ScrHight 0.85*ScrLength 0.75*ScrHight];

Color=['b' 'g' 'r' 'c' 'm' 'y' 'k' 'b' 'g' 'r' 'c' 'm' 'y' 'k']; % Same as MATLAB order
% Color=['r' 'b' 'g' 'k' 'm' 'c' 'y' 'r' 'b' 'g' 'k' 'm' 'c' 'y']; % Different from MATLAB order
% Shape=['o' 's' 'd' 'v' 'p' '*' 'x' '+' '.' '^' '<' '>' 'h'];
% [Color(i) Shape(i)])
MarkerSize=5;
set(gcf,'Position',figpos)

a = ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp}(:,cv_a);
b = ECVAm.Detail.ECVAmodel.ProjScores{no_of_comp}(:,cv_b);
CN = ECVAm.ClassmarkerNumeric;

plot(a,b,'w');

ix=zeros(max(CN),1);
for i=1:max(CN)
    ix_temp=find(CN==i);
    ix(i)=ix_temp(1);
end

if isempty(SampleNames) % Color according to class
    hold on
	for i=1:max(CN)
        h(i)=plot(a(CN==i),b(CN==i),'o','MarkerEdgeColor',Color(i),'MarkerFaceColor',Color(i),'MarkerSize',MarkerSize);
    end
    hold off
elseif strcmpi(SampleNames,'numbers') % Color according to class with numbers as labels
    number_string=int2str((1:length(CN))');
    hold on
	for i=1:max(CN)
        h(i)=plot(a(CN==i),b(CN==i),'o','MarkerEdgeColor',Color(i),'MarkerFaceColor',Color(i),'MarkerSize',MarkerSize);
        text(a(CN==i),b(CN==i),[ones(size(number_string(CN==i),1),1)*'  '  number_string(CN==i,:)],'Color',Color(i));
    end
    hold off
else % Color according to class with samplenames as labels
    hold on
	for i=1:max(CN)
        h(i)=plot(a(CN==i),b(CN==i),'o','MarkerEdgeColor',Color(i),'MarkerFaceColor',Color(i),'MarkerSize',MarkerSize);
        text(a(CN==i),b(CN==i),[ones(size(SampleNames(CN==i),1),1)*'  '  SampleNames(CN==i,:)],'Color',Color(i));
    end
    hold off
end

if isnumeric(ECVAm.ClassmarkerOriginal)
    legend(h,int2str(ECVAm.ClassmarkerOriginal(ix,:)));
else
    legend(h,ECVAm.ClassmarkerOriginal(ix,:));
end

sub_addcross;

s=sprintf('ECV# %g',cv_a);
xlabel(s)
s=sprintf('ECV# %g',cv_b);
ylabel(s)
title('Extended Canonical Variates')

function B = sub_simpls(X,Y,lv)

% SIMPLS function tailormade for ECVA
% NOTE: Y and t are NOT mean centered as they are in the original algorithm
%
% B = sub_simpls(X,Y,lv);

B = cell(lv,1);

% Y=mncn(Y);

S = X'*Y;
for i = 1:lv
  [~,~,q] = svd(S,0);
  q=q(:,1);
  r = S*q;
  t = X*r;
  % t = t-mean(t);
  normt = sqrt(t'*t);
  t = t/normt;
  r = r/normt;
  p = X'*t;
  q = Y'*t;
  % u = Y*q; % Not necessary when B is the only output
  v = p;
  if i > 1
    v = v - V*(V'*p);
    % u = u - T*(T'*u); % Not necessary when B is the only output
  end
  v = v/sqrt(v'*v);
  S = S - v*(v'*S);

  R(:,i) = r;
  % T(:,i) = t; % Not necessary when B is the only output
  % P(:,i) = p; % Not necessary when B is the only output
  Q(:,i) = q;
  % U(:,i) = u; % Not necessary when B is the only output
  V(:,i) = v;
end

for i = 1:lv
    B{i} = R(:,1:i)*Q(:,1:i)';
end

function Table = sub_classtable(ECVAm,no_of_comp)

% Input:
%   ECVAm: output from ecva('model')
%   no_of_comp [optional]: LVs in the PLS model. Default is the optimal found as minimum number of errors
% Output:
%   Table
% RB: Column 1 is sensitivity
% RB: Column 2 is specificity
% Table = sub_classtable(ECVAm,no_of_comp)

if nargin==1
    no_of_comp = ECVAm.OptComp;
end

UniqueGroups = unique(ECVAm.ClassmarkerNumeric);

ix=cell(length(UniqueGroups),1);
NoSamplesinGroups=zeros(length(UniqueGroups),1);
for i=1:length(UniqueGroups)
    ix{i}=find(ECVAm.ClassmarkerNumeric==i);
    NoSamplesinGroups(i) = size(ix{i},1);
end

InGroup=zeros(length(UniqueGroups)+1,length(UniqueGroups)+2); InGroup(1,1)=NaN;
InGroup(2:end,1)=UniqueGroups;
InGroup(2:end,2)=NoSamplesinGroups;
InGroup(1,3:end)=UniqueGroups; % Should be changed to ClassmarkerOriginal which may be a char
NotInGroup=zeros(length(UniqueGroups)+1,length(UniqueGroups)+2); NotInGroup(1,1)=NaN;
NotInGroup(2:end,1)=UniqueGroups;
NotInGroup(2:end,2)=NoSamplesinGroups;
NotInGroup(1,3:end)=UniqueGroups; % Should be changed to ClassmarkerOriginal which may be a char

for i=1:length(UniqueGroups)
    for j=1:length(UniqueGroups)
      InGroup(i+1,j+2) = sum(ECVAm.Detail.ClassPred{no_of_comp}(ix{i}) == UniqueGroups(j));
      NotInGroup(i+1,j+2) = sum(ECVAm.Detail.ClassPred{no_of_comp}(ix{i}) ~= UniqueGroups(j));
    end
end
Table.InGroupCell=num2cell(InGroup);
Table.InGroupCell{1,1}=' ';
Table.NotInGroupCell=num2cell(NotInGroup);

InGroupRed=InGroup;
InGroupRed(:,1:2)=[];
InGroupRed(1,:)=[];
TruePos=zeros(1,length(UniqueGroups));
for i=1:length(UniqueGroups)
     TruePos(i)=InGroupRed(i,i);
end
FalsePos=sum(InGroupRed) - TruePos;
FalseNeg=NoSamplesinGroups' - TruePos;
% TrueNeg=sum(NoSamplesinGroups) - TruePos - FalseNeg; % Old version
TrueNeg=sum(NoSamplesinGroups) - TruePos - FalseNeg - FalsePos; % Changed by RB

for i=1:length(UniqueGroups)
    Table.SensSpec(i,1) = 100*(TruePos(i)/(TruePos(i)+FalseNeg(i)));
    Table.SensSpec(i,2) = 100*(TrueNeg(i)/(TrueNeg(i)+FalsePos(i)));
end
Table.SensSpec = round(Table.SensSpec*10)/1000; % Three decimals
Table.ListOfMisClassifiedHeader=['SampleNo ' 'Original ' 'Estimated'];
Table.ListOfMisClassified=find(ECVAm.ClassmarkerNumeric~=ECVAm.Detail.ClassPred{no_of_comp});

ClassmarkerOriginalUnique = unique(ECVAm.ClassmarkerOriginal);
ClassmarkerNumericUnique = unique(ECVAm.ClassmarkerNumeric);

for i = 1:length(Table.ListOfMisClassified)
    ixtemp = ECVAm.Detail.ClassPred{no_of_comp}(Table.ListOfMisClassified(i)) == ClassmarkerNumericUnique;
    ClassmarkerOriginalPredicted(i) = ClassmarkerOriginalUnique(ixtemp);
end
if isempty(Table.ListOfMisClassified)
    Table.ListOfMisClassified = 'No misclassifications';
else
    Table.ListOfMisClassified=[Table.ListOfMisClassified ECVAm.ClassmarkerOriginal(Table.ListOfMisClassified) ClassmarkerOriginalPredicted'];
end

function sub_addcross

hold on
V=axis;
if sign(V(1))==sign(V(2)) && sign(V(3))==sign(V(4))
    % No lines through origo
elseif sign(V(1))==sign(V(2))
    plot([V(1) V(2)],[0 0],'Color',[0.7 0.7 0.7])
elseif sign(V(3))==sign(V(4))
    plot([0 0],[V(3) V(4)],'Color',[0.7 0.7 0.7])
else
    plot([V(1) V(2)],[0 0],[0 0],[V(3) V(4)],'Color',[0.7 0.7 0.7])
end
hold off

