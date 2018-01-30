function [w h distr]=lraNMF(v,opts)
%% v=w*h, where v, w ,h are nonnegative matrices.
%  usage:  >>opts=struct('NumOfComp',10);
%          >>[w h]=lraNMF(v,opts);
%   opts=struct('NumOfComp',[],'maxit',1000,'tol',1e-5,'alpha',0,'alg','mult','lra','pca','lra_rank',[],'w',[],'h',[],'trackit',20,'CNMF','none');
%        You can set some/all of them if necessary.
%   opts.NumOfComp: number of columns of w.
%       .alpha: parameter for sparseness/smoothness constraints
%       .alg:   mult|hals|apg   update formulas
%       .lra:   pca|random|randpca   low-rank approximation algorithms
%       .lra_rank: rank of low-rank approximations. lra_rank>=NumOfComp is
%                   recommended.
%       .w .h: initial guess of w and h.
%       .trackit: integer number for tracking the current status
%       .CNMF: a string referring to handle of a constrained NMF algorithm. h will be
%             refined by calling @CNMF(h).
%
% Please cite:
%    Guoxu Zhou; Cichocki, A.; Shengli Xie; , "Fast Nonnegative Matrix/Tensor Factorization Based on Low-Rank Approximation," 
%    IEEE Transactions on Signal Processing, vol.60, no.6, pp.2928-2940, June 2012
%    doi: 10.1109/TSP.2012.2190410
%    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166354&isnumber=6198804
%
%
%
% PLEASE FEEL FREE TO CONTACT ME FOR ANY QUESTIONS.
%
%                                by Guoxu Zhou
%                                zhoguoxu@brain.riken.jp
%                                http://bsp.brain.riken.jp
%                                Feb. 29, 2012.
%
%
% Refs for NMF algs:
% mult: D. D. Lee and H. S. Seung, ¡°Algorithms for non-negative matrix factorization,¡±
%       in Advances in Neural Information Processing Systems 13,
%       T. K. Leen, T. G. Dietterich, and V. Tresp, Eds. Cambridge, MA:
%       MIT Press, 2000, pp. 556¨C562
% HALS: A. Cichocki, R. Zdunek, A.-H. Phan, and S. Amari, Nonnegative Matrix
%     and Tensor Factorizations: Applications to ExploratoryMulti-Way
%     Data Analysis and Blind Source Separation. Chichester, U.K.:Wiley,
%     2009.
% APG [NeNMF]: Naiyang Guan; Dacheng Tao; Zhigang Luo; Bo Yuan; , "NeNMF: An Optimal Gradient 
%     Method for Nonnegative Matrix Factorization," IEEE Transactions on Signal Processing, 
%     vol.60, no.6, pp.2882-2898, June 2012

optsdef=struct('NumOfComp',[],'maxit',2000,'maxiniter',20,'tol',1e-5,'alpha',0,'alg','hals','lra','pca','lra_rank',[],'w',[],'h',[],'trackit',20,'CNMF','none');
if ~exist('opts','var')
    opts=struct;
end
[r,maxit,maxiniter,tol,alpha,mode,lra,lra_rank,w,h,trackit,CNMF]=scanparam(optsdef,opts);

if isempty(r)
    r=size(v,2);
end

if isempty(lra_rank)
    lra_rank=min(2*r,min(size(v)));
else
    lra_rank=min(max(lra_rank,r),min(size(v)));
end


if iscell(v)   % v=v{1}*v{2}
    if numel(v)~=2
        error('If v is a cell, only two elements are permitted.');
    else
        x=v{1};y=v{2};
        M=size(x,1);
        T=size(y,2);
    end
else
    [M T]=size(v);
    if M==r
        x=eye(M,M);
        y=v;
    else
        [x y]=lowrankapp(v,lra_rank,lra);
        if nargout>2
            distr=dist_fro(v,x,y);
        end
    end
end


if isempty(w) 
    w=rand(M,r); 
end
if isempty(h) 
    h=rand(r,T); 
end
switch lower(mode)
    case 'hals'
        maxit=maxit/r;
        cit=1;
        nh=sum(h.*h,2);
        for it=1:maxit
            w0=w;
            for i=1:r
                    w(:,i)=max(w(:,i)*nh(i)+x*(y*h(i,:)')-w*(h*h(i,:)'),eps)./(nh(i)+alpha);
                    nrm=max(sum(w(:,i).^2).^.5,eps);
                    w(:,i)=w(:,i)./nrm;
                    h(i,:)=nrm*h(i,:);

                    h(i,:)=max(h(i,:)+(w(:,i)'*x)*y-(w(:,i)'*w)*h,eps);
                    nh(i)=h(i,:)*h(i,:)';
            end
            
         if (trackit>0)&&(~rem(cit,trackit))
            if norm(w-w0,'fro')<tol
                break;
            end
        end           
           
        end
    case 'mult'
        for it=1:maxit
            w0=w;
            w=max(w.*(x*(y*h')./max(w*(h*h')+alpha.*w,eps)),eps);
%             w=w.*(max(x*(y*h'),eps)./max(w*(h*h')+alpha.*w,eps));
            wnorm=sum(w);
            w=bsxfun(@rdivide,w,max(wnorm,eps(wnorm)));
            gradH=(max((w'*x)*y,eps)./max((w'*w)*h,eps));
            h=max(h.*gradH,eps);     

%             fprintf('%d: res=%f - %f\n',it,norm(x*y-w*h,'fro'),norm(v-w*h,'fro'));
            if (trackit>0)&&(~rem(it,trackit))
                if norm(w-w0,'fro')<tol
                    break;
                end
            end
        end
    case 'apg'
        for it=1:maxit
            %% update w
            rho=1/norm(h'*h,'fro');
            alpha1=1;
            yk=w;
            hht=h*h';
            xyht=x*(y*h');
            for init=1:maxiniter
                w0=w;
                wg=-xyht+yk*hht;
                w=max(yk-rho*wg,1e-9);
                wdiff=w-w0;
                alpha0=alpha1;
                alpha1=(1+sqrt(4*alpha0^2+1))/2;
                yk=w+((alpha0-1)/alpha1)*wdiff;                
            end
            
            nrw=max(sum(w),1e-6);
            w=bsxfun(@rdivide,w,nrw);
            h=bsxfun(@times,h,nrw');
            
            %% update h
            rho=1/norm(w'*w,'fro');
            alpha1=1;
            wtw=w'*w;
            wtv=(w'*x)*y;
            yk=h;
            for init=1:maxiniter
                h0=h;
                hg=-wtv+wtw*yk;
                h=max(yk-rho*hg,1e-9);
                hdiff=h-h0;
                alpha0=alpha1;
                alpha1=(1+sqrt(4*alpha0^2+1))/2;
                yk=h+((alpha0-1)/alpha1)*hdiff;
            end
        end
    otherwise
        error('Unsuported parameter.');
end

if ~strcmp(CNMF,'none')
    if exist(CNMF,'file')
        cnmf=str2func(CNMF);
        [w2 h]=cnmf(h,r);
        w=w*w2;
    else
        fprintf('File [%s] does not exist.\n',CNMF);
    end
end
