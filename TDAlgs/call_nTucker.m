function [ Ycap,ExplX,Xm ] = call_nTucker( Y,opts )
%% Call Tucker decomposition in NWAY toolbox 3.2
% usage: Ycap=call_nTucker(Y,opts);
%
% opts.
% tol - Convergence criterion
%            The relative change in fit for which the algorithm stops.
%            Standard is 1e-6, but difficult data might require a lower value.
%
% initmode - Initialization method
%            This option is ignored if PARAFAC is started with old values.
%            If no default values are given the default Options(2) is 0.
%            The advantage of using DTLD or SVD for initialization is that
%            they often provide good starting values. However, since the
%            initial values are then fixed, repeating the fitting will give
%            the exact same solution. Therefore it is not possible to substantiate
%            if a local minimum has been reached. To avoid that use an initialization
%            based on random values (2).
%
%            0  = fit using DTLD/GRAM for initialization (default if
%                                 three-way and no missing and if sizes are
%                                 largere than number of factors at least
%                                 in two modes)
%            1  = fit using SVD vectors for initialization (default if higher than three-way or missing)
%            2  = fit using random orthogonalized values for initialization
%            10 = fit using the best-fitting models of several models
%            fitted using a few iterations
%
% Scaling
%            0 or 1 = default scaling (columns in mode one carry the variance)
%            2      = no scaling applied (hence fixed elements will not be modified
%
% showFit - How often to show fit
%            Determines how often the deviation between the model and the data
%            is shown. This is helpful for adjusting the output to the number
%            of iterations. Default is 10. If showfit is set to NaN, almost no
%            outputs are given
%
% maxiter - Maximal number of iterations
%            Maximal number of iterations allowed. Default is 2500.
%
% ConstrF  : Constraints that must apply to 'Factors'.
%            Define a row-vector of size N that describes how
%            each mode should be treated.
%            '0' orthogonality (default)
%            '1' non-negativity
%            '2' unconstrained
%            '4' unimodality and non-negativity.
%            E.g.: [0 2 1] yields ortho in first mode, uncon in the second
%            and non-neg in the third mode.
%            Note: The algorithm uses random values if there are no
%            non-negative components in the iteration intermediates. Thus,
%            if non-negativity is applied, the iterations may be
%            non-monotone in minor sequences.
% ConstrG  : Constraints that must apply to 'G'.
%            '[]' or '0' will not constrain the elements of 'G'.
%            To define what core elements should be allowed, give a core that
%            is 1 (one) on all active positions and zero elsewhere - this boolean
%            core array must have the same dimensions as defined by 'Fac'.
%
%TUCKER multi-way tucker model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: When estimating a Tucker model on data using non-orthogonal factors,
%       the sum of square of the core may differ between models of the
%       same dataset. This is in order since the factors may
%       thus be correlated. However, the expl. var. should always be the same.
%

% $ Version 2.01  $ May 2011 $ Fixed problem with Tucker2 calculations & missing 'Option'$ CA $ Not compiled $
% $ Version 2.003 $ Jan 2002 $ Fixed problem with length of factors under special conditions $ CA $ Not compiled $
% $ Version 2.002 $ Jan 2002 $ Fixed reshaping of old input G $ RB $ Not compiled $
% $ Version 2.001 $ July 2001 $ Changed problem with checking if Factors exist (should check if it exists in workspace specifically)$ RB $ Not compiled $
% $ Version 2.00  $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.12  $ Date 14. Nov. 1999 $ Not compiled $
%
%
% Copyright, 1998-, 2011-
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Claus A. Andersson
% Chemometrics Group - Life Sciences
% University of Copenhagen
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% E-mail: claus@andersson.dk // ca@life.ku.dk

defopts=struct('NumOfComp',[],'tol',1e-6,'initmode',0,'Scaling',1,'showFit',10,'maxiter',2500,...
    'constF',[],'constG','FileName','initialValues','FileName');
if ~exist('opts','var')
    opts=struct();
end
[Fac tol initmode Scaling showFit maxiter constF constGF initfile]=scanparam(defopts,opts);

% global TDALABHOME;
%     initf=horzcat(TDALABHOME,filesep,'userdata',filesep,strtrim(initfile));
initf=horzcat(strtrim(initfile));
fstatus=exist(initf,'file');
if fstatus~=2
    fprintf(2,'[TDALAB] File %s not found. Default values will be used instead.\n',constGF);
    G=[];
    Factors=[];
else
    W=load(initf);
    name=fieldnames(W);
    if strcmpi(class(W.(name{1})),'ttensor')
        G=double(W.(name{1}).core);
        Factors=W.(name{1}).U;
    else
        fprintf(2,'No valid variables found. Default values will be used instead.\n');
        G=[];
        Factors=[];
    end
end

constGF=horzcat(strtrim(constGF));
fstatus=exist(constGF,'file');
if fstatus~=2
    fprintf(2,'[TDALAB] File %s not found. Default values will be used instead.\n',constGF);
    constG=[];
else
    file=load(constGF);
    name=fieldnames(file);
    constG=file.(name{1});
end

opts=[tol initmode 0 Scaling showFit maxiter];

Y=double(Y);
if isempty(G)
    [Factors,G,ExplX,Xm]=tucker(Y,Fac,opts,constF,constG);
else
    [Factors,G,ExplX,Xm]=tucker(Y,Fac,opts,constF,constG,Factors,G);
end
Ycap=ttensor(tensor(G),Factors);
end

