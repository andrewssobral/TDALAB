function [ Ycap krpres] = mrcp( Y,opts )
%% N-way CP decomposition based on a proceeding 3-way decomposition
%  Usage: Ycap=n3cp(Y,opts);
%  opts
%      .NumOfComp: rank of Ycap
%      .modeGroup: group of the modes. Such as: 1 2 ; 3 ; 4 of a 4-way
%          tensor. Split by ';'.
%      .CP3_FILE: algorithm for 3-way CP
%      .lra: none|pca|sampling low-rank approximation method
%      .lra_rank: rank of lra (i.e., reduced dimension)
%      .krp_init: random|svd initialization method for Khatri-rao prod
%           approximation
%      .krp_maxit: maxiteration 
%      .krp_costol: stopping criterion 
%      .krp_constraints: none|nonnegative|sparse
%      .krp_cons_para: parameter which is used for krp_constraints.
defopts=struct('NumOfComp',[],'modeGroup','1;2;3','CP3_FILE','filename.mat','lra','pca','lra_rank',[],'krp_init','random','krp_maxit',50,'krp_costol',1e-2,...
    'krp_constraint','none','krp_cons_para',[],'dualDimR',false,'krpCheck',true);
if ~exist('opts','var')
    opts=struct();
end
[J, pmodestr,CP3_FILE,lra,lra_rank,krp_init,krp_maxit,krp_costol,krp_cons,krp_cons_para,dualDimR,krpCheck]=scanparam(defopts,opts);

Ydims=size(Y);
N=numel(Ydims);

if N<3
    error('Only works for high-way tensors.');
end
J=J(1);

%% parsing the pmodestr
pmodestr=regexp(pmodestr,';','split');
pmodestr=pmodestr(~cellfun(@(x) isempty(strtrim(x))', pmodestr));
pmode=cellfun(@str2num,pmodestr,'uni',false);
if any(cellfun('isempty',pmode))
    error('Matlab:ParamError','Wrong split of modes. Default settings are used instead.');
end
pmode{3}=setdiff(1:N,[pmode{1}, pmode{2}]);


%% partition the dims
Y3=permute(Y,[pmode{1} pmode{2} pmode{3}]);
Y3dims=[prod(Ydims(pmode{1})) prod(Ydims(pmode{2})) prod(Ydims(pmode{3}))];
Y3=tensor(reshape(Y3,Y3dims));
Y3init=Y3;

if prod(Y3dims(1:2))<J
    if ~strcmpi(lra,'none')
        warning('Rank deficient. lra=''none'' is used instead.');
        lra='none';
    end
elseif Y3dims(1)<J
    if ~strcmpi(lra,'none')
        if dualDimR
            warning('Rand deficient. dualDimR=false is used instead.');
            dualDimR=false;
        end
    end
end

%% construct a 3-way tensor
if dualDimR  %% using one-mode of Y3 to reconstruct the tensor
    switch lower(lra)
        case 'sampling'
            if isempty(lra_rank)
                lra_rank=repmat(5*J,1,2);
            elseif numel(lra_rank)==1
                lra_rank=repmat(lra_rank,1,2);
            end
            if lra_rank(1)<Y3dims(2)
                ord=randperm(Y3dims(2));
                sel=sort(ord(1:lra_rank(1)));
                Y3=Y3(:,sel,:);
            end
            if lra_rank(2)<Y3dims(3)
                ord=randperm(Y3dims(3));
                sel=sort(ord(1:lra_rank(2)));
                Y3=Y3(:,:,sel);
            end
        case 'pca'
            if isempty(lra_rank)
                lra_rank=repmat(J,1,2);
            elseif numel(lra_rank)==1
                lra_rank=repmat(lra_rank,1,2);
            end
            lra_rank=max(lra_rank,J);
            lra_rank=min(lra_rank,Y3dims(2:3));
            U=cell(1,2);
            U{2}=nvecs(Y3,3,lra_rank(2));
            
            for it=1:5
                G=ttm(Y3,U{2},3,'t');
                U{1}=nvecs(G,2,lra_rank(1));
                G=ttm(Y3,U{1},2,'t');
                U{2}=nvecs(G,3,lra_rank(2));
            end
            Y3=ttm(Y3,U,[2 3],'t');
    end
else   %% dualDimR=False
    switch lower(lra)
        case 'sampling'
            if isempty(lra_rank)
                lra_rank=5*J;
            end
            if lra_rank<Y3dims(3)
                ord=randperm(Y3dims(3));
                sel=sort(ord(1:lra_rank));
                Y3=Y3(:,:,sel);
            end
        case 'pca'
            if isempty(lra_rank)
                lra_rank=J;
            end
            if size(Y3,3)>2*J
%                 U=nvecs(Y3,3,lra_rank(end));
                
                temp=double(tenmat(Y3,3));
                U=lowrankapp(temp,lra_rank(end),'pca');
                clear temp;
                
                Y3=ttm(Y3,U,3,'t');
            end
    end
end


%% low-order CP
% % for TDALAB
% global TDALABHOME;
% CP3_FILE=horzcat(TDALABHOME,filesep,'userdata',filesep,strtrim(CP3_FILE));

fstatus=exist(CP3_FILE,'file');
if fstatus~=2
    fprintf(2,'[TDALAB] File %s does not exist.\n',CP3_FILE);
    error('The algorithm is terminated automatically.');
end
y3alg=load(CP3_FILE);
if iscell(y3alg.algname)
    y3alg.algname=y3alg.algname{1};
    y3alg.alg=y3alg.alg{1};
    y3alg.algopts=y3alg.algopts{1};
end
% if isfield(y3alg.algopts,'NumOfComp')&&J~=y3alg.algopts.NumOfComp
%     fprintf(2,'[CPN3]: in %s the NumOfComp is %d (should be %d). \n',y3alg.algname,y3alg.algopts.NumOfComp(1),J);
    y3alg.algopts.NumOfComp=J;
% end

Y3=y3alg.alg(Y3,y3alg.algopts);
Y3.U{3}=bsxfun(@times,Y3.U{3},Y3.lambda');
Y3.lambda=ones(J,1);

minKrpApp=1;

krp_opts=struct('init',krp_init,'powit',krp_maxit,'constraints',krp_cons,'cons_para',krp_cons_para,'costol',krp_costol);
% n1cp_opts=struct('krp_init',krp_init,'krp_maxit',krp_maxit,'krp_costol',krp_costol,'constraints',krp_cons,'cons_para',krp_cons_para);
P=cell(1,N);    
if strcmpi(lra,'none')  %% without low-rank approximation
    for idx=1:3
        if numel(pmode{idx})>1
            [P(pmode{idx})]=pikrpapprox(Y3.U{idx},Ydims(pmode{idx}),krp_opts);
        else
            P{pmode{idx}}=Y3.U{idx};
        end
    end
else    
    if dualDimR
        %% pmode{1}
        if numel(pmode{1})>1
            [P(pmode{1})]=pikrpapprox(Y3.U{1},Ydims(pmode{1}),krp_opts);
            if krpCheck
                krpApp=1-norm(Y3.U{1}-khatrirao(P(pmode{1}),'r'),'fro')/norm(Y3.U{1},'fro');
                minKrpApp=min(krpApp,minKrpApp);
            end
        else
            P{pmode{1}}=Y3.U{1};
        end
        %% resmode
        Yp=ttm(Y3init,(Y3.U{1}'*Y3.U{1})\Y3.U{1}',1);
        Yp=double(tenmat(Yp,1))';
%         Yp=(Y3.U{1}\double(tenmat(Y3init,1)))';
        resmode=[pmode{2} pmode{3}];
        [P(resmode)]=pikrpapprox(Yp,Ydims(resmode),krp_opts);
        if krpCheck
            krpApp=1-norm(Yp-khatrirao(P(resmode),'r'),'fro')/norm(Yp,'fro');
            minKrpApp=min(krpApp,minKrpApp);
        end
    else
        for idx=1:2
            if numel(pmode{idx})>1
                [P(pmode{idx})]=pikrpapprox(Y3.U{idx},Ydims(pmode{idx}),krp_opts);
                if krpCheck
                    krpApp=1-norm(Y3.U{idx}-khatrirao(P(pmode{idx}),'r'),'fro')/norm(Y3.U{idx},'fro');
                    minKrpApp=min(minKrpApp,krpApp);
                end
            else
                P{pmode{idx}}=Y3.U{idx};
            end 
        end
            %% pmode 3
        U=khatrirao(Y3.U{2},Y3.U{1})';
        Yp=double(tenmat(Y3init,3));
        Yp=Yp/U;
        [P(pmode{3})]=pikrpapprox(Yp,Ydims(pmode{3}),krp_opts);  
        if krpCheck
            krpApp=1-norm(Yp-khatrirao(P(pmode{3}),'r'),'fro')/norm(Yp,'fro');
            minKrpApp=min(minKrpApp,krpApp);
        end
    end
end
Ycap=ktensor(P);    
if krpCheck
    fprintf(2,'[MRCPD] Worst Khatri-Rao Fit (KRP Fit) is %f\n',minKrpApp);
    if minKrpApp<.5
        fprintf(2,'[MRCPD] Bad KRP Fit. Please change the settings or try another algorithm.\n');
    elseif minKrpApp<.9
        fprintf(2,'[MRCPD] The KRP Fit looks not so good. Perhaps change the settings or try another algorithm.\n');
    end
end
end

