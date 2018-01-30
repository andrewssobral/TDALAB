function [w h]=semiNMF(v,r,opts)
%% coded by Guoxu Zhou
% zhouguoxu@brain.riken.jp
optsdef=struct('maxit',1000,'tol',1e-5,'alpha',0,'alg','semi','lra','pca','w0',[],'h0',[]);
if ~exist('opts','var')
    opts=struct;
end
[maxit,tol,alpha,mode,lra,w,h]=scanparam(optsdef,opts);



[M T]=size(v);
if M==r
    x=eye(M,M);
    y=v;
else
    [x y]=lowrankapp(v,r,lra);
end

if isempty(w) 
    w=randn(M,r); 
end
if isempty(h) 
    h=rand(r,T); 
end

switch mode
    case 'cov'
    case 'semi'
        cit=1;
        for it=1:maxit
            w0=w;
            w=x*(y/h);
                        
            XtF=(w'*x)*y;
            FtF=w'*w;
            XtFp=max(XtF,0);
            XtFm=XtFp-XtF;
            FtFp=max(FtF,0);
            FtFm=FtFp-FtF;  
            h=h.*(max(XtFp+FtFm*h,eps)./max(XtFm+FtFp*h,eps)).^0.5;
            
%             if norm(w0-w,'fro')<tol
%                 break;
%             end

        end
    otherwise
        error('Unsuported parameter.');
end
