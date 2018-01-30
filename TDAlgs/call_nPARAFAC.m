function Ycap=call_nPARAFAC(X,opts)
%% This function calls parafac.m included in N-WAY toolbox 3.2
%   Usage: Ycap=call_nPARAFAC(X,opts);
%
%opts.
% tol - Convergence criterion
%            The relative change in fit for which the algorithm stops.
%            Standard is 1e-6, but difficult data might require a lower value.
%
% initmode - Initialization method
%            This option is ignored if PARAFAC is started with old values.
%            If no default values are given the default value is 0.
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
% const      A vector telling type of constraints put on the loadings of the
%            different modes. Same size as DimX but the i'th element tells
%            what constraint is on that mode.
%            0 => no constraint,
%            1 => orthogonality
%            2 => nonnegativity
%            3 => unimodality (and nonnegativitiy)
%            4 => L1 fitting (will be imposed in all modes)
%            5 => L1 fitting and nonnegativity (will be imposed in all modes)
%            If const is not defined, no constraints are used.
%            For no constraints in a threeway problem const = [0 0 0]
% FixMode    FixMode is a binary vector of same sixe as DimX. If
%            FixMode(i) = 1 => Mode i is fixed (requires old values given)
%            FixMode(i) = 0 => Mode i is not fixed hence estimated
%            Ex.: FixMode = [0 1 1] find the scores of a data set given the loadings.
%            When some modes are fixed, the numbering of the components will
%            also be fixed. Normally components are sorted according to variance
%            as in PCA, but this will not be performed if some modes are fixed.
%
% Weights    If a matrix of the same size as X is given, weighted regression
%            is performed using the weights in the matrix Weights. Statistically
%            the weights will usually contain the inverse error standard
%            deviation of the particular element
%      >>>>>> Please specify a mat-file which defines a matrix named as
%            'Weights'
%
%
%     C. A. Andersson and R. Bro
%     The N-way Toolbox for MATLAB
%     Chemometrics & Intelligent Laboratory Systems. 52 (1):1-4, 2000.
%     http://www.models.life.ku.dk/source/nwaytoolbox/
%
% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA  02110-1301, USA.



% $ Version 1.03 $ Date 1. October   1998 $ Not compiled $ Changed sign-convention because of problems with centered data
% $ Version 1.04 $ Date 18. February 1999 $ Not compiled $ Removed auxiliary line
% $ Version 1.06 $ Date 1. December  1999 $ Not compiled $ Fixed bug in low fit error handling
% $ Version 1.07 $ Date 17. January  2000 $ Not compiled $ Fixed bug in nnls handling so that the algorithm is not stopped until nonnegative appear
% $ Version 1.08 $ Date 21. January  2000 $ Not compiled $ Changed init DTLD so that primarily negative loadings are reflected if possible
% $ Version 1.09 $ Date 30. May 2000 $ Not compiled $ changed name noptioPF to noptiopf
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.001 $ June 2001 $ Fixed error in weighted regression $ RB $ Not compiled $
% $ Version 2.002 $ Jan 2002 $ Fixed scaling problem due to non-identifiability of DTLD(QZ) by scaling and normalizing after each iteration $ RB $ Not compiled $
% $ Version 2.003 $ Jan 2002 $ Fixed negative solutions when nonneg imposed $ RB $ Not compiled $
% $ Version 2.004 $ Jan 2002 $ Changed initialization when many components used $ RB $ Not compiled $
% $ Version 2.005 $ Jan 2002 $ Changed absolute fit criterion (approacing eps) into relative sse/ssx$ RB $ Not compiled $
% $ Version 2.006 $ Jan 2002 $ Fixed post-scaling when fixed loadings $ RB $ Not compiled $
% $ Version 2.01 $ Jan 2003 $ Removed corcondia for two-way data (doesn't work) and fixed a bug for data with dimension 2 $ RB $ Not compiled $
% $ Version 2.011 $ feb 2003 $ Added an option (4) for not post scaling components $ RB $ Not compiled $
% $ Version 2.10  $ jan 2004 $ Fixed a plotting error occuring when fitting model to old data $ RB $ Not compiled $
% $ Version 2.11  $ jan 2004 $ Fixed that PCA can be fitted $ RB $ Not compiled $
% $ Version 2.12  $ Jul 2004 $ Fixed initialization bug $ RB $ Not compiled $
% $ Version 2.13 $ Jan 2005 $ Modified sign conventions of scores and loads $ RB $ Not compiled $
% $ Version 2.14 $ Feb 2006 $ Fixed bug in sign-swicth when loadings are fixed $ RB $ Not compiled $
% $ Version 2.15 $ Aug 2007 $ ALS scheme substantially improved using linesearch acc. to Rajih, Comon, Harshman 2007 $ RB $ Not compiled $
% $ Version 2.16 $ Jun 2010 $ Fixed that scaling of component is not done in fixed modes $ RB $ Not compiled $


defopts=struct('NumOfComp',[],'tol',1e-6,'initmode',0,'Scaling',0,'showFit',NaN,'maxiter',2500,...
    'const',[],'FixMode',[],'Weights','FileName','OldLoad','FileName');

if ~exist('opts','var')
    opts=struct();
end

[FAC tol initmode Scaling showFit maxiter const  FixMode weifile oldLoad]=scanparam(defopts,opts);


fstatus=exist(weifile,'file');
if fstatus~=2
    fprintf(2,'[TDALAB] File %s for ''Weights'' not found. Default values will be used instead.\n',weifile);
    Weights=[];
else
    W=load(weifile);
    fns=fieldnames(W);
    Weights=W.(fns{1});    
end

fstatus=exist(oldLoad,'file');
if fstatus~=2
    fprintf(2,'[TDALAB] File %s for ''OldLoad'' not found. Default value will be used instead.\n',oldLoad);
    oldLoad=[];
else
    oldLoad=load(oldLoad);
    ns=fieldnames(oldLoad);
    if iscell(oldLoad.(ns{1}))
        oldLoad=oldLoad.(ns{1});
    elseif strcmpi(class(oldLoad.(ns{1})),'ktensor')
        oldLoad=oldLoad.(ns{1}).U;
    else
        fprintf(2,'No valid variables found. Default value will be used instead.\n');
        oldLoad=[];
    end
end


parafac_opts=[tol initmode 0 Scaling showFit maxiter];

X=double(X);
FAC=FAC(1);
Ycap=parafac(X,FAC,parafac_opts,const,oldLoad,FixMode,Weights);
Ycap=ktensor(ones(FAC,1),Ycap);
