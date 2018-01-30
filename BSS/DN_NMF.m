function [A,X,res] = DN_NMF(Y,opts)
%function [A,X,res] = nmf_qp_cg(Y,A,X,MaxIter,tol_c,lambda)
% The Damped Newton (DN) algorithm developed by
% Copyrirght R. Zdunek, A.-H. Phan, and A. Cichocki
%
% DN

%% line 61 is changed
defopts = struct('NumOfComp',[],'A0',[],'X0',[],'MaxIter',500,'Tol',1e-6,'lambda',1);  % Theta maps checkStep.
if ~exist('opts','var')
    opts = struct;
end

[r,A,X,MaxIter,tol_c,lambda] ...
    =scanparam(defopts,opts);
if isempty(r)
    r=size(Y,2);
end

[I,T] = size(Y);
epsil = eps;
Y = max(eps,Y);
n = 0; k = 0;

if isempty(A)
    A=rand(I,r);
end
if isempty(X)
    X=rand(r,T);
end

% Number of components
r = size(A,2);

% Scaling
da = sqrt(sum(A.^2,1))+eps;
A = bsxfun(@rdivide,A,da);
X = bsxfun(@times,da',X);

    k = 0;
    res = [];    c = [];
    while k < MaxIter

        k = k + 1;
%         if ~mod(k,50)
%             disp(['Iteration: ',num2str(k)]);
%         end
      
% UPDATES FOR A 
% ======================================================================
        % Memory prelocation    
        M = I*r;
        Q = spalloc(M,M,M*r); Q_tilde = spalloc(M,M,M*r);
        B = X*X'; C = X*Y';
        Q = kron(speye(I),B); % Hessian of D(Y||AX) with respect to A
        Qd = zeros(M,1);
        At = A';   at = At(:); atp = at;
        GAt = B*A' - C;

        active = find((at <= epsil) & (GAt(:) > 0)); % active set
        degradated = find((at <= epsil) & (abs(GAt(:)) < epsil), 1);
        if ~isempty(degradated)
           disp(['Degenerate: ',num2str(k)]);
        end
        inactive = setdiff(1:length(at),active); % inactive set

        if ~isempty(active)
            a_tilde = at;            a_tilde(active) = 0;
            c_tilde = C(:);          c_tilde(active) = 0;
            Q_tilde = Q;
            Q_tilde(:,active) = 0;   Q_tilde(active,:) = 0;
            Qd(active) = 1;
            Q_tilde = Q_tilde + spdiags(Qd,0,M,M) + lambda*speye(M);
        else
            a_tilde = at;
            Q_tilde = Q + lambda*speye(M);     
            c_tilde = C(:);
        end

        % CGS 
        tol = 1e-6 + exp(-k);
        at = zeros(M,1); Io = [];
        Ao = setdiff(1:M,Io);   Qd = zeros(M,1);
        
        while ~isempty(Ao)
              [a_tilde,FLAG,RELRES,ITER] = pcg(Q_tilde,c_tilde,tol);
              Ao = find(a_tilde < 0);  
              Q_tilde(Ao,:) = 0; Q_tilde(:,Ao) = 0; Q_tilde = Q_tilde + spdiags(Qd,0,M,M);
              c_tilde(Ao) = 0; 
        end
        at = a_tilde; 

        At = (max(epsil,reshape(at,r,I)));
        A = At';
        
        da = eps+sqrt(sum(A.^2,1));
        A = bsxfun(@rdivide,A,da);
        X = bsxfun(@times,da',X);
        
% UPDATES FOR X 
% ======================================================================
        
        HX = A'*A + 1e-12*eye(r);
        GX = HX*X - A'*Y;
        I_active = find((X < eps)&(GX > 0));
        UX = GX;   
        UX = HX\UX;  UX(I_active) = 0; 
        X = max(eps,X - UX);
        
        % Normalization
%         da = eps+sqrt(sum(A.^2,1));
%         A = bsxfun(@rdivide,A,da);
%         X = bsxfun(@times,da',X);
%         
        dx = eps+sqrt(sum(X.^2,2));
        X = bsxfun(@ldivide,dx,X);
        A = bsxfun(@times,dx',A);

        res(k) = norm(Y - A*X,'fro')/norm(Y,'fro');

 % Stagnation breaking       
        if k > 3
            c(k) = res(k-1) - res(k);
            if (c(k)) < tol_c
                lambda = lambda/2;
            end
            if k > 30 & c(k) < 0 & c(k-1) < 0 & c(k-2) < 0
                break
            end
        end
                    
end % while k
     




