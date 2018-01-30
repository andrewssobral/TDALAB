function [ Q Y] = lowrankapp( X, r,alg,p)
%% min ||X-Q*Q'*X|| s.t. Q'*Q=I_r.
%        r<<M,N
% default_p=2;

%% PCA is considered first
if strcmpi(alg,'pca')
    [Q d V]=svds(X,r,'L');
    flag=diag(d)>1e-6;
    Q=Q(:,flag);
    Y=Q'*X;
    return;
end

%% other algorithms
N=size(X,2);
if nargin<=3
    p=default_p;
end
if isempty(p)
    p=default_p;
end


switch lower(alg)
    case {'sampling','rs'}
        X0=X;
        X=X(max(abs(X))>eps,:);
        N=size(X,2);
        L=min(p*r,N);
        o=randperm(N);
        X=X(:,o(1:L));
        if r<L
            [Q d]=lowrankapp(X,r,'pca');
            Y=Q'*X0;
        else
            Q=X;
            Y=Q\X0;
        end
    case {'randpca','random'}
        L=max(p,1)*r;
        if 2*L<N
            Y=X*randn(N,L);        
            Q=lowrankapp(Y,r,'pca');
        else
            Q=lowrankapp(X,r,'pca');
        end
        Y=Q'*X;
    otherwise
        error('Unsuported algorithm.');
end





% err=norm(X-Q*Y,'fro');
% disp(time);
% disp(err);

