function [w h]=rslraNMF(v,opts)
optsdef=struct('NumOfComp',[],'maxit',1000,'tol',1e-5,'alpha',0,'alg','mult','lra','pca','trackit',20,'SamplingPara',1);
if ~exist('opts','var')
    opts=struct;
end
[r,maxit,tol,alpha,mode,lra,trackit,rs]=scanparam(optsdef,opts);

[M T]=size(v);


if isempty(r)
    r=M;
end
if M==r
    x=eye(M,M);
    y=v;
elseif T==r
    x=v;
    y=eye(T,T);
else
    [x y]=lowrankapp(v,r,lra);
end


SAPR=0.9;
rind=1:M;cind=1:T;
if rs>1
    p=rs*r;
    if p<T*SAPR
        cind=randperm(T);    
        cind=cind(1:min(T,p));
    end
    if p<M*SAPR
        rind=randperm(M);
        rind=rind(1:min(M,p));
    end
elseif (rs>0)&&(rs<SAPR)
    rmax=max(r,ceil(rs*M));
    cmax=max(r,ceil(rs*T));
    cind=randperm(T);
    cind=cind(1:cmax);
    rind=randperm(M);
    rind=rind(1:rmax);
end


    w=rand(M,r); 
    h=rand(r,T); 

mode=lower(mode);    
switch mode
    case 'hals'
        maxit=maxit/r;
        cit=1;
        for it=1:maxit
            w0=w;
            for i=1:r
                left=[1:i-1 i+1:r];
                h(i,:)=max((w(rind,i)'*x(rind,:)*y-w(rind,i)'*w(rind,left)*h(left,:)),eps)/(w(rind,i)'*w(rind,i));
                w(:,i)=max((x*(y(:,cind)*h(i,cind)')-w(:,left)*(h(left,cind)*h(i,cind)')),eps)/(h(i,:)*h(i,:)'+alpha);

                cit=cit+1;
                if (cit>maxit*0.2)&&(~rem(cit,trackit))
                    if norm(w-w0,'fro')<tol
                        break;
                    end
                end
            end
            nw=sum(w);
            w=bsxfun(@rdivide,w,max(nw,eps(nw)));
            h=bsxfun(@times,h,nw');
        end
    case 'mult'
        for it=1:maxit
            w0=w;h0=h;
            w=w.*(max(x*(y(:,cind)*h(:,cind)'),eps)./max(w*(h(:,cind)*h(:,cind)')+alpha.*w,eps));
            wnorm=sum(w);
            w=bsxfun(@rdivide,w,max(wnorm,eps(wnorm)));
            gradH=(max(w(rind,:)'*x(rind,:)*y,eps)./max(w(rind,:)'*w(rind,:)*h,eps));
            h=h.*gradH;
            
            if (it>maxit*0.2)&&(~rem(it,trackit))
                if norm(w-w0,'fro')<tol
                    break;
                end
            end
        end
    otherwise
        error('Unsuported parameter.');
end
% 
% if ~strcmp(CNMF,'none')
%     if exist(CNMF,'file')
%         cnmf=str2func(CNMF);
%         [w2 h]=cnmf(h,r);
%         w=w*w2;
%     else
%         fprintf('File [%s] does not exist.\n',CNMF);
%     end
% end
