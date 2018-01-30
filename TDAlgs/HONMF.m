function [Core, FACT, varexpl] = HONMF(X, d, varargin)
% HONMF: Higher Order NON-NEGATIVE MATRIX FACTORIZATION
%
% Authors:
%   Morten Mørup
%   Technical University of Denmark,
%   Institute for Matematical Modelling
%
% Reference:
%   M. Mørup et al. Algorithms for Sparse Higher Order Non-negative Tensor Factorization. Technical University of Denmark, 2006.
%
% Usage:
%       [Core, FACT] = HONMF(X, d, varargin)
%
% Example:
%       opts.lambda=[1 0 0 1];
%       opts.maxiter=1000;
%       X=rand(10,20,30);
%       [Core, FACT] = HONMF(X, [3 3 3], opts)
%
% Input:
%   X                 Array of non negative data
%   d                 1 x length(FACT) vector containing number of factors
%                     for each modality, i.e. d(i) is number of factors in
%                     FACT{i}.
%   options
%     .costfcn        Cost function to optimize
%                       'ls': Least squares (default)
%                       'kl': Kullback Leibler
%     .FACT           Initial FACT
%     .Core           Initial Core
%     .W              Indicator tensor of same size as X, 0 where data is missing
%                     1 where present
%     .constFact      1 x length(FACT) vector indicating which factors to update,
%                     i.e. constFACT(i)=0: updates FACT{i}  (default)
%                          constFACT(i)=1: FACT{i} is not updated.
%     .constCore      constCore=0: updates Core  (default)
%                     constCore=1: Core is not updated.
%     .lambda         Sparsity weight on Core and Factors,
%                     (sparsity function used is the 1-norm)
%                       lambda(1):      Core
%                       lambda(i+1):    FACT{i}
%     .maxiter        Maximum number of iterations (default 100)
%     .conv_criteria  Function exits when cost/delta_cost exceeds this
%     .accel          Wild driver accelleration parameter (default 1.3)
%     .displaylevel   Level of display: [off | final | iter]
%
% Output:
%   Core                 d(1) x d(2) x ... x d(N) array
%   FACT                 cell array of factors for each modality, i.e.
%                        FACT{i} is a size(V,i) x d(i) matrix
%   varexpl              Explained variation
%
% Copyright (C) Morten Mørup and Technical University of Denmark,
% September 2006
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%   Edit history:
%   30 Oct 2008:
%           varexplained changed from relative to
%           SST=norm(X-mean(X),'fro')^2 to SST= norm(X,'fro')^2
%           error corrected in initialization

% -------------------------------------------------------------------------
% Parse input arguments
if nargin>=3, opts = varargin{1}; else opts = struct; end
costfcn = mgetopt(opts, 'costfcn', 'ls', 'instrset', {'ls','kl'});
constFACT=mgetopt(opts,'constFACT',zeros(size(d)));
constCore=mgetopt(opts,'constCore',0);
lambda = mgetopt(opts, 'lambda', zeros(1,1+ndims(X)));
mX=max(X(:));
for i=1:ndims(X)
    FACT{i}=mX*rand(size(X,i),d(i));
    if lambda(i+1)==0 && sum(lambda)~=0
        FACT{i}=normalize(FACT{i});
    end
end
FACT = mgetopt(opts, 'FACT', FACT);
Core = mgetopt(opts, 'Core', rand(d));
if sum(lambda(2:end))>0 && lambda(1)==0
    Core=Core/sqrt(sum(Core(:).^2)+eps);
end
nuF = mgetopt(opts, 'nuF', 1);
nuC = mgetopt(opts, 'nuC', 1);
W = mgetopt(opts,'W',ones(size(X)));
X=W.*X;
maxiter = mgetopt(opts, 'maxiter', 2500);
conv_criteria = mgetopt(opts, 'conv_criteria', 1e-6);
accel = mgetopt(opts, 'accel', 1.3);
beta = mgetopt(opts, 'beta', 2);
displaylevel = mgetopt(opts, 'displaylevel', 'iter', 'instrset', ...
    {'off','iter','final'});

% -------------------------------------------------------------------------
% Initialization
sst = sum(X(:).^2);
Rec = honmf_rec(Core,FACT);
sse = 2*costLS(X,W.*Rec);
sparse_cost=0;

for i=1:length(FACT)
    sparse_cost=sparse_cost+lambda(i+1)*sum(FACT{i}(:));
end
switch costfcn
    case 'ls'
        cost= .5*sse + sparse_cost+lambda(1)*sum(Core(:));
    case 'kl'
        xlogxx=X.*log(X+eps)-X;
        xlogx_x=sum(xlogxx(:));
        cost =  costKL(X,W.*Rec,xlogx_x)+sparse_cost+ lambda(1)*sum(Core(:));
end
delta_cost = 1;
iter = 0;
keepgoing = 1;
cost_old=cost;


% -------------------------------------------------------------------------
% Display information
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s','Iteration','Expl. var.','nuC','nuF','Cost func.','Delta costf.',' Time(S)');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+--------------+');
if any(strcmp(displaylevel, {'final','iter'}))
    disp('Higher Order Non-negative Matrix Factorization');
    disp(['To stop algorithm press ctrl-C'])
    disp(['Sparsity set to ' num2str(lambda)])
    disp(['Constant Core ' num2str(constCore) ' Constant Fact ' num2str(constFACT)]);
    if sum(W(:))==numel(W)
        disp(['No values missing'])
    else
        disp(['Values of X being zero treated as missing, a total of : ' num2str(length(W(W==0))/prod(size(W))*100) ' %' ])
    end
    disp(['Optimization method ' costfcn])
    disp('');

end

t=cputime;
% -------------------------------------------------------------------------
% Optimization loop
while keepgoing
    pause(0.001);
    told=t;
    if mod(iter,100)==0
        if any(strcmp(displaylevel, {'iter'}))
            disp(dline); disp(dheader); disp(dline);
        end
    end

    % Update Core and FACT
    switch costfcn
        case 'ls'
            [Core, FACT, cost, sse, Rec, nuF, nuC] = ...
                sls_update(X,Core, FACT, cost, sse, Rec,constCore,constFACT,lambda,nuF,nuC,accel,beta,W);
        case 'kl'
            [Core, FACT, cost, sse, Rec, nuF, nuC] = ...
                skl_update(X,Core, FACT, cost, sse, Rec,constCore,constFACT,lambda,nuF,nuC,accel,beta,xlogx_x,W);
    end

    delta_cost = cost_old - cost;
    cost_old=cost;
    iter=iter+1;
    t=cputime;

    % Display information
    %    if any(strcmp(displaylevel, {'iter'}))
    if mod(iter,5)==0
        disp(sprintf('%12.0f | %12.4f | %12.6f | %12.6f | %12.6f |  %12.6f |  %12.6f ', ...
            iter,(sst-sse)/sst,nuC,nuF,cost,delta_cost,t-told));
    end

    % Check if we should stop
    if delta_cost<cost*conv_criteria || cost<10e-10*sst
        % Small improvement with small step-size
        if nuF<=accel && nuC<=accel
            if any(strcmp(displaylevel, {'iter','final'}))
                disp('HONMF has converged');
            end
            keepgoing = 0;
            % Small improvement - maybe because of too large step-size?
        else
            nuF = 1;
            nuC = 1;
        end
    end
    % Reached maximum number of iterations
    if iter>=maxiter
        if any(strcmp(displaylevel, {'iter','final'}))
            disp('Maximum number of iterations reached');
        end
        keepgoing=0;
    end
end
varexpl=(sst-sse)/sst;

% -------------------------------------------------------------------------
% Sparse Least squares update function
function [Core, FACT, cost, sse, Rec, nuF, nuC] = ...
    sls_update(X,Core, FACT, cost, sse, Rec,constCore,constFACT,lambda,nuF,nuC,accel,beta,W);

pause(0.00001); %Enables to break algorithm by pushing "Control C".
g=find(lambda>0);
le=length(FACT);
if constCore==0   % Update Core
    cont=0;
    cost_old=cost;
    Coreold=Core;
    Recold=Rec;
    while cont==0
        sparse_Cost=0;
        Dx=X;
        Dy=W.*Rec;
        for i=1:length(FACT)
            sparse_Cost=sparse_Cost+lambda(i+1)*sum(FACT{i}(:));
            Dx=tmult(Dx,FACT{i}',i);
            Dy=tmult(Dy,FACT{i}',i);
        end
        if ~isempty(g) && lambda(1)==0 && length(g)~=length(lambda)  % If some of the modalities are sparsified but not the current, use normalized cost function. If all modalities are sparsified, don't normalize core.
            tx = Core.*Dy;
            tx =sum(tx(:));
            ty = Core.*Dx;
            ty=sum(ty(:));
            Dx = Dx + Core*tx;
            Dy = Dy + Core*ty;
        end
        Core=Core.*(((Dx+eps)./(Dy+ lambda(1)+eps)).^nuC);
        if ~isempty(g) && lambda(1)==0 && length(g)~=length(lambda)
            Core=Core/sqrt(sum(Core(:).^2)+eps);
        end
        Rec=reconstruct(Core,FACT);
        sse=2*costLS(X,W.*Rec);
        cost=0.5*sse+lambda(1)*sum(Core(:))+sparse_Cost;
        if cost<=cost_old
            nuC=accel*nuC;
            cont=1;
        else
            nuC=max([nuC/beta,1]);
            cont=0;
            Core=Coreold;
            Rec=Recold;
        end
    end
end

cont=0;
cost_old=cost;
FACTold=FACT;
Recold=Rec;
while cont==0
    sparse_cost=0;
    for i=1:le
        if constFACT(i)==0 % Update FACT
            ind=1:le;
            ind=setdiff(ind,i);
            Z=Core;
            for k=ind
                Z=tmult(Z,FACT{k},k);
            end
            Z=matrizicing(Z,i);
            Xi=matrizicing(X,i);
            Dx=Xi*Z';
            Wi=matrizicing(W,i);
            Dy=(Wi.*(FACT{i}*Z))*Z';
            if lambda(i+1)==0 && ~isempty(g)   % If some of the modalities are sparsified but not the current, use normalized cost function
                tx = sum(Dy.*FACT{i},1);
                ty = sum(Dx.*FACT{i},1);
                Dx = Dx + repmat(tx,[size(FACT{i},1),1]).*FACT{i};
                Dy = Dy + repmat(ty,[size(FACT{i},1),1]).*FACT{i};
            end
            FACT{i}=FACT{i}.*(((Dx+eps)./(Dy+eps+lambda(i+1))).^nuF);
            if lambda(i+1)==0 && ~isempty(g)
                FACT{i}=normalize(FACT{i});
            end
            Rec=unmatrizicing(FACT{i}*Z,i,size(X));
        end
        sparse_cost=sparse_cost+lambda(i+1)*sum(FACT{i}(:));
    end
    sse=2*costLS(X,W.*Rec);
    cost=0.5*sse+lambda(1)*sum(Core(:))+sparse_cost;
    if cost<=cost_old
        nuF=nuF*accel;
        cont=1;
    else
        nuF=max([nuF/beta,1]);
        cont=0;
        FACT=FACTold;
        Rec=Recold;
    end
end

if sum(lambda)==0 % No sparseness
    [Core, FACT]=scale(Core,FACT);
end





% -------------------------------------------------------------------------
% Sparse Kullback Leibler update function
function [Core, FACT, cost, sse, Rec, nuF, nuC] = ...
    skl_update(X,Core, FACT, cost, sse, Rec,constCore,constFACT,lambda,nuF,nuC,accel,beta,xlogx_x,W);

pause(0.00001); %Enables to break algorithm by pushing "Control C".
g=find(lambda>0);
le=length(FACT);
if constCore==0   % Update Core
    cont=0;
    cost_old=cost;
    Coreold=Core;
    Recold=Rec;
    while cont==0
        sparse_Cost=0;
        Dx=X./(Rec+eps);
        Dy=W;
        for i=1:length(FACT)
            sparse_Cost=sparse_Cost+lambda(i+1)*sum(FACT{i}(:));
            Dx=tmult(Dx,FACT{i}',i);
            Dy=tmult(Dy,FACT{i}',i);
        end
        if ~isempty(g) && lambda(1)==0 && length(g)~=length(lambda)  % If some of the modalities are sparsified but not the current, use normalized cost function. If all modalities are sparsified, don't normalize core.
            tx = Core.*Dy;
            tx =sum(tx(:));
            ty = Core.*Dx;
            ty=sum(ty(:));
            Dx = Dx + Core*tx;
            Dy = Dy + Core*ty;
        end
        Core=Core.*(((Dx+eps)./(Dy+ lambda(1)+eps)).^nuC);
        if ~isempty(g) && lambda(1)==0 && length(g)~=length(lambda)
            Core=Core/sqrt(sum(Core(:).^2)+eps);
        end
        Rec=reconstruct(Core,FACT);
        cost=costKL(X,W.*Rec,xlogx_x)+lambda(1)*sum(Core(:))+sparse_Cost;
        if cost<=cost_old
            nuC=accel*nuC;
            cont=1;
        else
            nuC=max([nuC/beta,1]);
            cont=0;
            Core=Coreold;
            Rec=Recold;
        end
    end
end
cont=0;
cost_old=cost;
FACTold=FACT;
Recold=Rec;
while cont==0
    sparse_cost=0;
    for i=1:le
        if constFACT(i)==0 % Update FACT
            ind=1:le;
            ind=setdiff(ind,i);
            Z=Core;
            for k=ind
                Z=tmult(Z,FACT{k},k);
            end
            Z=matrizicing(Z,i);
            Xi=matrizicing(X,i);
            Reci=matrizicing(Rec,i);
            Dx=(Xi./(Reci+eps))*Z';
            Dy=matrizicing(W,i)*Z';
            if lambda(i+1)==0 && ~isempty(g) % If some of the modalities are sparsified but not the current, use normalized cost function
                tx = sum(Dy.*FACT{i},1);
                ty = sum(Dx.*FACT{i},1);
                Dx = Dx + repmat(tx,[size(FACT{i},1),1]).*FACT{i};
                Dy = Dy + repmat(ty,[size(FACT{i},1),1]).*FACT{i};
            end
            FACT{i}=FACT{i}.*(((Dx+eps)./(Dy+eps+lambda(i+1))).^nuF);
            if lambda(i+1)==0 && ~isempty(g)
                FACT{i}=normalize(FACT{i});
            end
            sparse_cost=sparse_cost+lambda(i+1)*sum(FACT{i}(:));
            Rec=unmatrizicing(FACT{i}*Z,i,size(X));
        end
    end
    cost=costKL(X,W.*Rec,xlogx_x)+lambda(1)*sum(Core(:))+sparse_cost;
    if cost<=cost_old
        nuF=nuF*accel;
        cont=1;
    else
        nuF=max([nuF/beta,1]);
        cont=0;
        FACT=FACTold;
        Rec=Recold;
    end
end

if sum(lambda)==0 % No sparseness
    [Core, FACT]=scale(Core,FACT);
end
sse=2*costLS(X,W.*Rec);


% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname);
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end


% -------------------------------------------------------------------------
% Normalize W
function F = normalize(F)
Q = sqrt(sum(F.^2,1));
F = F./repmat(Q+eps,[size(F,1),1]);


% -------------------------------------------------------------------------
% Scale Core and FACT to ensure W and H keep within reasonable value range
function [Core, FACT] = scale(Core, FACT)
for i=1:length(FACT)
    Q = sqrt(sum(FACT{i}.^2,1));
    FACT{i} = FACT{i}.*repmat(1./(Q+eps),[size(FACT{i},1),1]);
    Core=tmult(Core,diag(Q),i);
end


% -------------------------------------------------------------------------
% Reconstructs the data from Core and FACT
function Rec = reconstruct(Core, FACT)
Rec=Core;
for i=1:length(FACT)
    Rec=tmult(Rec,FACT{i},i);
end



% -------------------------------------------------------------------------
function cost=costLS(V,Rec,missing)
cost =0.5*norm(V(:)-Rec(:),'fro')^2;


% -------------------------------------------------------------------------
function cost=costKL(V,Rec,xlogx_x,W,missing)
cost = (-V.*log(Rec+eps)+Rec);
cost=sum(cost(:))+xlogx_x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Other functions
function Rec = honmf_rec(Core,FACT)
Rec=Core;
for i=1:length(FACT)
    Rec=tmult(Rec,FACT{i},i);
end

function A=krprod(B,C)
sb=size(B,1);
sc=size(C,1);
A=zeros(sb*sc,size(B,2));
for k=1:size(B,2)
    A(:,k)=reshape(C(:,k)*B(:,k)',sb*sc,1);
end

function Y=matrizicing(X,n)
N=ndims(X);
Y=reshape(permute(X, [n 1:n-1 n+1:N]),size(X,n),numel(X)/size(X,n));

function T=outerprod(FACT)
T=0;
for i=1:size(FACT{1},2)
    Y=1;
    for j=1:length(FACT)
        U=FACT{j};
        Y=tmult(Y,U(:,i),j);
    end
    T=T+Y;
end

function A=tmult(T,M,n)
Dt=size(T);
Dm=size(M);
Tn=matrizicing(T,n);
Tnew=M*Tn;
Dt(n)=Dm(1);
A=unmatrizicing(Tnew,n,Dt);

function X=unmatrizicing(X,n,D)
if n==1
    perm=[1:length(D)];
else
    perm=[2:n 1 n+1:length(D)];
end

X=permute(reshape(X,D([n 1:n-1 n+1:length(D)])),perm);