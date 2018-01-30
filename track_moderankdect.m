function [rk rkt sap eigt err]=track_moderankdect(X,n,SamplingPara,Er,truerank)
%% mode rank detection based on SORTE using randomly sampled fibers
X=double(X);
Xdim=size(X);
NumOfMode=numel(Xdim);
DIM=prod(Xdim([1:n-1 n+1:end]));
        
        if DIM<500
            z=double(tenmat(tensor(X),n));
            rk=sorte(z);
            return;
        end
 p0=10;           
        NumOfSampled=10;
       
        z=zeros(NumOfSampled,Xdim(n));
        index=1:NumOfMode;index(n)=[];
        var=cell(1,NumOfMode);
        var{n}=1:Xdim(n);
 track_it=1;      
        o=randperm(DIM);
        curri=0;rk0=0;
        for t=1:DIM
            [var{index}]=ind2sub(Xdim([1:n-1 n+1:end]),o(t));
            zt=reshape(double(X(var{:})),1,Xdim(n));
            if norm(zt)>0.001
                curri=curri+1;
                z(curri,:)=zt;
            end
fprintf('curri=%d *** \n',curri);
            if ~rem(curri,SamplingPara)
                [eigvect eigv]=eig(cov(z(1:curri,:)'));
                [eigv ind]=sort(diag(eigv),'descend');
                rk=sorte(eigv);
                curri
                rkt(track_it)=rk;
                sap(track_it)=curri;
                eigt{track_it}=eigv;
                if curri<=truerank
                    err(track_it)=(norm(X-z'*(z'\X),'fro')-Er)/norm(X,'fro');
                else
                    [uz dz]=eigs(cov(z'),truerank);
                    uz=uz'*z;
                    err(track_it)=(norm(X-uz'*(uz'\X),'fro')-Er)/norm(X,'fro');
                end
%                 err(track_it)=norm(X-z'*inv(z*z')*z*X,'fro')/Er;
                track_it=track_it+1;
                
                if (curri>Xdim(n)-2*SamplingPara)
                    break;
                else
                    rk0=rk;
                end
            end
        end
        fprintf('# rank=%d   curr=%d\n',rk,curri);
        
end