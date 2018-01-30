function [ P ] = n1cp( Y,pmode,U_pmode,opts )
Ydim=size(Y);
N=numel(Ydim);
defopts=struct('krp_init','random','krp_maxit',50,'krp_costol',1e-2,'constraints','none','cons_para',0);
if ~exist('opts','var')
    opts=struct();
end
[krp_init,krp_maxit,krp_costol,cons,cons_para]=scanparam(defopts,opts);

if strcmp(class(Y),'tensor')~=1
    Y=tensor(Y);
end

if isempty(pmode)
    error('pmode is not specified');
end

P=cell(N,1);
np=numel(pmode);

if iscell(U_pmode)
    J=size(U_pmode{1},2);
    if numel(U_pmode)~=np
        error('Size of pmode and U_pmode does not match.');
    end    
    P(pmode)=U_pmode;
else
    J=size(U_pmode,2);
    if np~=1
        error('Size of pmode and U_pmode does not match.');
    end
    P{pmode}=U_pmode;
end

% %% update Y
if np>1
    alldims=1:N;
    resdims=setdiff(alldims,pmode);
    Ipmode=prod(Ydim(pmode));
    Ires=prod(Ydim(resdims));  

    %% way - I --
    rkra=khatrirao(U_pmode,'r');
    Yp=double(tenmat(Y,resdims));
    Yp=Yp/rkra';
else
    resdims=[1:pmode-1 pmode+1:N];
    Yp=ttm(Y,(U_pmode'*U_pmode)\U_pmode',pmode);
    Yp=tenmat(Yp,pmode);
    Yp=double(Yp)';
end

%% update Y


%% Set up and error checking on initial guess for U.

%% reset the dimension
resYdims=Ydim(resdims);
if numel(resYdims)>1
    krp_opts=struct('init',krp_init,'powit',krp_maxit,'constraints',cons,'cons_para',cons_para,'costol',krp_costol);
    [P(resdims)]=pikrpapprox(Yp,resYdims,krp_opts);
else
    P{resdims}=Yp;
end

%% improve



%% return the results
P=ktensor(ones(J,1),P);

end