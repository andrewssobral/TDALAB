function algsInitialize()
%% Important!!!
fprintf('[TDALAB] Loading and filtering algorithms ... ');
global algs paraList paraTypeList;

ind=1;
algs(ind).name='call_ASD';algs(ind).details='ASD algorithm';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','Tol',1e-10,'lambda',1e-3,'MaxIter',10000);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'Tol',paraType.paraDouble,'lambda',paraType.paraDouble,'MaxIter',paraType.paraDouble);

ind=ind+1;
algs(ind).name='swatld';algs(ind).details='swatld algorithm';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','MaxIter',500,'Tol',1e-9);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_APTLD';algs(ind).details='APTLD algorithm';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','MaxIter',500,'Tol',1e-9,'p',1e+20,'q',1e+20,'r',1e+20);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'p',paraType.paraDouble,'q',paraType.paraDouble,'r',paraType.paraDouble);

%% CP-OMP
ind=ind+1;
algs(ind).name='cp_smbss';algs(ind).details='[CP-SMBSS] CP based on One Single Mode BSS';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','BSSmode',1,'BSSalgFile','Click here to choose','lra','pca|sampling|randpca|none','lra_rank',[],'krp_init','random|svd',...
    'krp_maxit',50,'krp_costol',1e-2,'constraints','none|nonnegative','cons_para',0,'verbose',true);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'BSSmode',paraType.paraDouble,'BSSalgFile',paraType.paraFile,'lra',paraType.paraMenu,'lra_rank',paraType.paraDouble,'krp_init',paraType.paraMenu,...
    'krp_maxit',paraType.paraDouble,'krp_costol',paraType.paraDouble,'constraints',paraType.paraMenu,'cons_para',paraType.paraDouble,'verbose',paraType.paraTF);


%% CP-OMP
ind=ind+1;
algs(ind).name='mrcp';algs(ind).details='[MRCPD] CP Decomposition using Mode Reduction';
algs(ind).nonnegativity=AlgType.Both;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','modeGroup','1 ; 2 ; 3 4','CP3_FILE','fileName.mat','lra','pca|sampling|none','lra_rank',[],'krp_init','random|svd','krp_maxit',50,'krp_costol',1e-2,...
    'krp_constraint','none|nonnegative|sparse','krp_cons_para',0,'dualDimR',false,'krpCheck',true);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'modeGroup',paraType.paraString,'CP3_FILE',paraType.paraFile,'lra',paraType.paraMenu,'lra_rank',paraType.paraDouble,'krp_init',paraType.paraMenu,'krp_maxit',paraType.paraDouble,'krp_costol',paraType.paraDouble,...
    'krp_constraint',paraType.paraMenu,'krp_cons_para',paraType.paraDouble,'dualDimR',paraType.paraTF,'krpCheck',paraType.paraTF);


ind=ind+1;
algs(ind).name='lraHALSCP';algs(ind).details='HALS CP Decomposition';
algs(ind).nonnegativity=AlgType.Both;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Rank of output','lra_rank',[],'lra_iter',50,'lra_model','CP|Tucker|none','initmode','CP|random|sampling','maxiter',50,'inner_maxiter',20,'tol',1e-6,'constraint','none|nonnegative|sparse','cons_param',[]);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'lra_rank',paraType.paraDouble,'lra_iter',paraType.paraDouble,...
    'lra_model',paraType.paraMenu,'initmode',paraType.paraMenu,'maxiter',paraType.paraDouble,'inner_maxiter',paraType.paraDouble,'tol',paraType.paraDouble,'constraint',paraType.paraMenu,'cons_param',paraType.paraDouble);



ind=ind+1;
algs(ind).name='lroat';algs(ind).details='LROAT: Low Rank Orthogonal Approximation of Tensors.';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','MaxIter',100,'Tol',1e-4,'init','nvecs|random|load');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'init',paraType.paraMenu);


ind=ind+1;
algs(ind).name='call_cp_als';algs(ind).details='CP_ALS (tensor toolbox 2.5)';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','tol',1e-12,'maxiters',500,...
    'init','random|nvecs|load','printitn',10,'dimorder','1:N');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'tol',paraType.paraDouble,'maxiters',paraType.paraDouble,...
    'init',paraType.paraMenu,'printitn',paraType.paraDouble,'dimorder',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_tucker_als';algs(ind).details='Tucker ALS [HOOI] (tensor toolbox 2.5)';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','tol',1e-12,'maxiters',50,...
    'init','random|nvecs|load','printitn',10,'dimorder','1:N');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'tol',paraType.paraDouble,'maxiters',paraType.paraDouble,...
    'init',paraType.paraMenu,'printitn',paraType.paraDouble,'dimorder',paraType.paraDouble);



ind=ind+1;
algs(ind).name='call_cp_nmu';algs(ind).details='Nonnegative CP_ALS (tensor toolbox 2.5)';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','tol',1e-12,'maxiters',500,...
    'init','random|nvecs|load','printitn',10,'dimorder','1:N');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'tol',paraType.paraDouble,'maxiters',paraType.paraDouble,...
    'init',paraType.paraMenu,'printitn',paraType.paraDouble,'dimorder',paraType.paraDouble);


ind=ind+1;
algs(ind).name='FSTDT1opt';algs(ind).details='FSTD1 algorithm';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('Rank','rank for sampling','MaxIter',1000,'epsilon',1e-9,'Verbose',false,'checkStep',10);
paraTypeList{ind}=struct('Rank',paraType.paraDouble,'MaxIter',paraType.paraDouble,'epsilon',paraType.paraDouble,'Verbose',paraType.paraTF,'checkStep',paraType.paraDouble);


ind=ind+1;
algs(ind).name='FSTDT2rnd';algs(ind).details='FSTD2 algorithm';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('Rank','rank for sampling');
paraTypeList{ind}=struct('Rank',paraType.paraDouble);


ind=ind+1;
algs(ind).name='FastNTFAPG';algs(ind).details='<html><font color=''blue''>[FastNTF_APG] Fast NTF Using Accelerated Proximal Gradient</f></html>';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=Inf;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Rank of output','lra_rank','rank for low-rank approximation','CPalg','FileName','FastMode',true,'initmode','CP|random','maxiter',500,'max_initer',20,'tol',1e-6,...
    'alpha',1,'beta',0.005);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'lra_rank',paraType.paraDouble,'CPalg',paraType.paraFile,'FastMode',paraType.paraTF,'initmode',paraType.paraMenu,...
    'maxiter',paraType.paraDouble,'max_initer',paraType.paraDouble,'tol',paraType.paraDouble,...
    'alpha',paraType.paraDouble,'beta',paraType.paraDouble);


ind=ind+1;
algs(ind).name='call_HONMF';algs(ind).details='High Order NMF (HONMF)';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Yes;
paraList{ind}=struct('NumOfComp','Columns of each mode','maxiter',100,'conv_criteria',1e-6,...
    'constFACT','zero vector','constCore',0,'lambda','Sparsity weight on core and factors','accel',1.3,'costfcn','ls|kl','displaylevel','off|final|iter');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'maxiter',paraType.paraDouble,'conv_criteria',paraType.paraDouble,...
    'constFACT',paraType.paraDouble,'constCore',paraType.paraDouble,'lambda',paraType.paraDouble,'accel',paraType.paraDouble,'costfcn',paraType.paraMenu,'displaylevel',paraType.paraMenu);


ind=ind+1;
algs(ind).name='cp_alsls';algs(ind).details='CP_ALSLS for 3,4,5-order tensor';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=5;
algs(ind).model=AlgType.CP;
algs(ind).acceptTensor=AlgType.Yes;
paraList{ind}=struct('NumOfComp',0,'Tol',[1e-6 1e-4],'MaxIter',[50 5000],'lsearch','lsb|lsh|none|elsr|elsc',...
    'Ninit',3,'comp','off|on','Ainit','valid only for 3-way tensor|load');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'Tol',paraType.paraDouble,'MaxIter',paraType.paraDouble,'lsearch',paraType.paraMenu,...
    'Ninit',paraType.paraDouble,'comp',paraType.paraMenu,'Ainit',paraType.paraMenu);


ind=ind+1;
algs(ind).name='appr3d';algs(ind).details='Fast Tucker approximation [appr3d]';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Yes;
paraList{ind}=struct('MaxIter',500,'Tol',1e-9,'Verbose',false);
paraTypeList{ind}=struct('MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'Verbose',paraType.paraTF);


ind=ind+1;
algs(ind).name='BSSTucker';algs(ind).details='<html><font color=''blue''>MBSS_Tucker_ALS (MBSS)</font></html>';
algs(ind).nonnegativity=AlgType.Both;
algs(ind).numOfMode=Inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Yes;
paraList{ind}=struct('NumOfComp',0,'MaxIter',50,'DimOrder',[],'norm',2,...
    'PMFalgIDs','PMF algorithm ID for each mode','SamplingPara',10,'SamplingThres',1e-2,'SamplingMode','times|FSTD1|FSTD2','Verbose',false,'update',true);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'DimOrder',paraType.paraDouble,'norm',paraType.paraDouble,...
    'PMFalgIDs',paraType.paraPMF,'SamplingPara',paraType.paraDouble,'SamplingThres',paraType.paraDouble,'SamplingMode',paraType.paraMenu,'Verbose',paraType.paraTF,'update',paraType.paraTF);



ind=ind+1;
algs(ind).name='stmlsvd';algs(ind).details='Sequentially truncated multilinear SVD (ST-MSVD)';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=Inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Columns of each mode','p','Order of modes');
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'p',paraType.paraDouble);


ind=ind+1;
algs(ind).name='call_bcdLMN';algs(ind).details='Block component decomposition: bcdLMN_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('R',[],'L',[],'M',[],'N',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('R',paraType.paraDouble,'L',paraType.paraDouble,'M',paraType.paraDouble,'N',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_bcdLrMrNr';algs(ind).details='Block component decomposition: bcdLrMrNr_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('L_vec',[],'M_vec',[],'N_vec',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('L_vec',paraType.paraDouble,'M_vec',paraType.paraDouble,'N_vec',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_bcdLL1';algs(ind).details='Block component decomposition: bcdLL1_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('R',[],'L',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('R',paraType.paraDouble,'L',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_bcdLM';algs(ind).details='Block component decomposition: bcdLM_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('R',[],'L',[],'M',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('R',paraType.paraDouble,'L',paraType.paraDouble,'M',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_bcdLrLr1';algs(ind).details='Block component decomposition: bcdLrLr1_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('L_vec',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('L_vec',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);

ind=ind+1;
algs(ind).name='call_bcdLrMr';algs(ind).details='Block component decomposition: bcdLrMr_alsls';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=3;
algs(ind).model=AlgType.BCD;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('L_vec',[],'M_vec',[],'lsearch','elsr|none|lsh|lsb|elsc','comp','on|off',...
    'Tol1',1e-6,'MaxIt1',500,'Tol2',1e-4,'MaxIt2',50,'Ninit',3);
paraTypeList{ind}=struct('L_vec',paraType.paraDouble,'M_vec',paraType.paraDouble,'lsearch',paraType.paraMenu,'comp',paraType.paraMenu,...
    'Tol1',paraType.paraDouble,'MaxIt1',paraType.paraDouble,'Tol2',paraType.paraDouble,'MaxIt2',paraType.paraDouble,'Ninit',paraType.paraDouble);




ind=ind+1;
algs(ind).name='lraSNTD';algs(ind).details='<html><font color=''blue''>[lraSNTD]Sequential Nonnegative Tucker Decomposition based on lraNMF</f></html>';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=Inf;
algs(ind).model=AlgType.Tucker;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'maxiter',5000,'maxiniter',10,'dimorder',[],'tol',1e-6,'trackit',20,'norm',1,'alpha',0,'SamplingPara',1,...
    'lra','pca|random|randpca','alg','HALS|apg|mult','update',true,'nncore',true);
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'maxiter',paraType.paraDouble,...
    'maxiniter',paraType.paraDouble,'dimorder',paraType.paraDouble,'tol',paraType.paraDouble,'trackit',paraType.paraDouble,'norm',paraType.paraDouble,'alpha',paraType.paraDouble,'SamplingPara',paraType.paraDouble,...
    'lra',paraType.paraMenu,'alg',paraType.paraMenu,'update',paraType.paraTF,'nncore',paraType.paraTF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% functions for n-way toolbox
if exist('parafac.m','file')
    fs=which('parafac.m');
    fs=regexp(fs,filesep,'split');
    if strcmpi('nway',fs{end-1}(1:4))
        %% add the functions of nway
        
        %% parafac
        ind=ind+1;
        algs(ind).name='call_nPARAFAC';algs(ind).details='PARAFAC (N-Way toolbox 3.20)';
        algs(ind).nonnegativity=AlgType.Both;
        algs(ind).numOfMode=Inf;
        algs(ind).model=AlgType.CP;
        algs(ind).acceptTensor=AlgType.Both;
        paraList{ind}=struct('NumOfComp','Rand of output CP tensor','tol',1e-6,'initmode',0,'Scaling',0,'showFit',10,'maxiter',2500,...
            'const',0,'FixMode',0,'Weights','FileName.mat','OldLoad','FileName.mat(cell)');
        paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'tol',paraType.paraDouble,'initmode',paraType.paraDouble,'Scaling',paraType.paraDouble,'showFit',paraType.paraDouble,'maxiter',paraType.paraDouble,...
            'const',paraType.paraDouble,'FixMode',paraType.paraDouble,'Weights',paraType.paraFile,'OldLoad',paraType.paraFile);

        %% Tucker
        ind=ind+1;
        algs(ind).name='call_nTucker';algs(ind).details='Tucker (N-Way toolbox 3.20)';
        algs(ind).nonnegativity=AlgType.Both;
        algs(ind).numOfMode=Inf;
        algs(ind).model=AlgType.Tucker;
        algs(ind).acceptTensor=AlgType.Both;
        paraList{ind}=struct('NumOfComp','FAC','tol',1e-6,'initmode',0,'Scaling',1,'showFit',10,'maxiter',2500,...
    'constF',[],'constG','FileName','initialValues','FileName');
        paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'tol',paraType.paraDouble,'initmode',paraType.paraDouble,'Scaling',paraType.paraDouble,'showFit',paraType.paraDouble,'maxiter',paraType.paraDouble,...
    'constF',paraType.paraDouble,'constG',paraType.paraFile,'initialValues',paraType.paraFile);

        
    end
end
%% end of N-way   


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ===== PMF ===== %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Efficient FastICA
algs(ind).name='PMFefica';algs(ind).details='Efficient FastICA';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Number of Components/sources','initmode','randn|eye','SaddleTest',true);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'initmode',paraType.paraMenu,'SaddleTest',paraType.paraTF);

%% Efficient Robust ICA
ind=ind+1;
algs(ind).name='PMFerica';algs(ind).details='Equivariant Robust ICA';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Number of Components/sources','MaxIter',1000,'Tol',1e-3,'whitening',false);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'whitening',paraType.paraTF);

%%
%% SOBI algorithm
ind=ind+1;
algs(ind).name='PMFsobi';algs(ind).details='SOBI: Second-Order Blind Identification';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Number of Components/sources','p','number of correlation matrices');
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'p',paraType.paraDouble);

%% WASOBI algorithm
ind=ind+1;
algs(ind).name='PMFwasobi';algs(ind).details='WASOBI: Weights-Adjusted SOBI';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Number of Components/sources','AR_order',1,'rmax',1);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'AR_order',paraType.paraDouble,'rmax',paraType.paraDouble);

%% nmfsc algorithm
ind=ind+1;
algs(ind).name='call_ICA_EBM';algs(ind).details='ICA by Entropy Bound Minimization';
algs(ind).nonnegativity=AlgType.No;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'lra','none|pca|sampling|random','lra_p',[],'max_iter_fastica',100,'max_iter_orth',1000,...
    'max_iter_orth_refine',1000,'max_iter_nonorth',1000,'max_cost_increase_number',5,...
    'stochastic_search_factor',1,'saddle_test_enable',true,'verbose',true);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'lra',paraType.paraMenu,'lra_p',paraType.paraDouble,'max_iter_fastica',paraType.paraDouble,'max_iter_orth',paraType.paraDouble,...
    'max_iter_orth_refine',paraType.paraDouble,'max_iter_nonorth',paraType.paraDouble,'max_cost_increase_number',paraType.paraDouble,...
    'stochastic_search_factor',paraType.paraDouble,'saddle_test_enable',paraType.paraTF,'verbose',paraType.paraTF);


%% lraNMF algorithm
ind=ind+1;
algs(ind).name='lraNMF';algs(ind).details='[lraNMF] NMF based on low-rank approximation';
algs(ind).nonnegativity=AlgType.Both;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'maxit',5000,'maxiniter',20,'tol',1e-5,'alpha',0,'trackit',20,'alg','hals|apg|mult','lra','pca|random|randpca','lra_rank',[]);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'maxit',paraType.paraDouble,...
    'maxiniter',paraType.paraDouble,'tol',paraType.paraDouble,'alpha',paraType.paraDouble,'trackit',paraType.paraDouble,'alg',paraType.paraMenu,'lra',paraType.paraMenu,'lra_rank',paraType.paraDouble);


%% VCA-NMF
ind=ind+1;
algs(ind).name='nmfVCA';algs(ind).details='NMF based on Vertex Comp. Analysis (VCA)';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'SNR',[],'verbose',false);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'SNR',paraType.paraDouble,'verbose',paraType.paraTF);

%% DNNMF algorithm
ind=ind+1;
algs(ind).name='DN_NMF';algs(ind).details='Damped Newton (DN) algorithm';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp','Number of Components/sources','MaxIter',500,'Tol',1e-5,'lambda',1);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'MaxIter',paraType.paraDouble,'Tol',paraType.paraDouble,'lambda',paraType.paraDouble);

%% NNMF algorithm in MATLAB
ind=ind+1;
algs(ind).name='PMFnnmf';algs(ind).details='NNMF (provided by MATLAB)';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'algorithm','als|mult','Display','off|final|iter','MaxIter',100,'TolFun',1e-4,'TolX',1e-4,...
    'UseParallel','never|always','UseSubstreams','never|always','replicates',1);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'algorithm',paraType.paraMenu,'Display',paraType.paraMenu,'MaxIter',paraType.paraDouble,'TolFun',paraType.paraDouble,'TolX',paraType.paraDouble,...
    'UseParallel',paraType.paraMenu,'UseSubstreams',paraType.paraMenu,'replicates',paraType.paraDouble);


%% beta_NMF_ME algorithm
ind=ind+1;
algs(ind).name='beta_nmf_ME';algs(ind).details='ME algorithm for NMF with the beta-divergence';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'beta',0.5,'n_iter',200);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'beta',paraType.paraDouble,'n_iter',paraType.paraDouble);

%% smoISNMF algorithm
ind=ind+1;
algs(ind).name='PMFSmoISNMF';algs(ind).details='Smooth Itakura-Saito NMF algorithm';
algs(ind).nonnegativity=AlgType.Yes;
algs(ind).numOfMode=2;
algs(ind).model=AlgType.PMF;
algs(ind).acceptTensor=AlgType.Both;
paraList{ind}=struct('NumOfComp',[],'n_iter',200,'lambda',10);
%%-----------------------------------------------------------------------------
paraTypeList{ind}=struct('NumOfComp',paraType.paraDouble,'n_iter',paraType.paraDouble,'lambda',paraType.paraDouble);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Add a new algorithm after this line.
    

%% Add your algorithm BEFORE this line.
%% ----------------------------------------------------
%% END OF THIS FILE. DON'T CHANGE ANYTHING AFTER HERE.
fprintf('Complete. \n');
end

