function [Y]=appr3d(X,opts)
%Fast Tucker approximation
% Approximation is sought in form of a = g x u1 x u2 x u3,
% Tucker decomposition
%
% Ref:
% Tucker Dimensionality Reduction of Three-Dimensional Arrays in Linear Time
% SIAM. J. Matrix Anal. & Appl. Volume 30, Issue 3, pp. 939-956 (2008) 
% Published September 25, 2008 
% 
defopts = struct('Tol',1e-10,'MaxIter',200,'Verbose',false);
if ~exist('opts','var')
    opts = struct;
end
[eps,niters,verbose]=scanparam(defopts,opts);
fprintf('Warning! The rank of loading matrix is detected automatically!\n');
fprintf('Please wait for processing now ...\n');
X=tensor(X);
sizeX=size(X);
n1=sizeX(1);
n2=sizeX(2);
n3=sizeX(3);
mat.n1=sqrt(n1);mat.n2=sqrt(n2);mat.n3=sqrt(n3);
[u1,u2,u3,G]=appr3d(@getele31,mat,n1,n2,n3,eps,niters,verbose);
G=tensor(G);
Y=ttensor(G,u1,u2,u3);
fprintf('The result dimension of the tucker core: %d x %d x %d\n',size(G));


    function [el]=getele31(i,j,k,mat)
        el=X(i,j,k);
    end


    function [u1,u2,u3,core]=appr3d(f,mat,n1,n2,n3,eps,niters,verbose)
        r1=2;
        r2=2;
        r3=2;
        u1=random('norm',0,1,[n1,r1]);
        u2=random('norm',0,1,[n2,r2]);
        u3=random('norm',0,1,[n3,r3]);
        %%
        %fm=zeros(n1,n2,n3);
        %for i1=1:n1
        %    for i2=1:n2
        %        for i3=1:n3
        %            fm(i1,i2,i3)=f(i1,i2,i3,mat);
        %        end
        %    end
        %end
        %%
        nrm=0;
        for iter=1:niters
            if ( iter ~= 1 )
                u1=[u1 uadd1];
                u2=[u2 uadd2];
                u3=[u3 uadd3];
            end
            
            [u1,tmp]=qr(u1,0);
            [u2,tmp]=qr(u2,0);
            [u3,tmp]=qr(u3,0);
            [n1,r1]=size(u1);
            [n2,r2]=size(u2);
            [n3,r3]=size(u3);
            ind1=1:r1;
            ind2=1:r2;
            ind3=1:r3;
            ind1=maxvol2(u1,ind1);
            ind2=maxvol2(u2,ind2);
            ind3=maxvol2(u3,ind3);
            %Compute core: A \approx U P V^t
            core=zeros(r1,r2,r3);
            for s1=1:r1
                for s2=1:r2
                    for s3=1:r3
                        core(s1,s2,s3)=f(ind1(s1),ind2(s2),ind3(s3),mat);
                    end
                end
            end
            uu1=inv(u1(ind1,:)');
            uu2=inv(u2(ind2,:)');
            uu3=inv(u3(ind3,:)');
            core=convt(core,uu1,uu2,uu3);
            oldnorm=nrm;
            nrm=norm(core(:));
            %Compress submatrix
            
            [u01,u02,u03,core]=svd3l(core,eps);
            nrm1=norm(core(:));
            u1=u1*u01;
            u2=u2*u02;
            u3=u3*u03;
            [r1 r2 r3]=size(core);
            %%
            %appr=convt(core,u1',u2',u3');
            %appr1=appr-fm;
            %fprintf('true residue: %e \n',norm(appr1(:))/norm(fm(:)));
            %%
            if verbose~=0
                fprintf('norm=%f,dnrm=%e,ranks=%d %d %d\n',nrm,abs(nrm-oldnorm)/abs(nrm),r1,r2,r3);
            end
            %Add columns to u and v
            [n1,r1]=size(u1);
            [n2,r2]=size(u2);
            [n3,r3]=size(u3);
            uadd1=zeros(n1,r1);
            uadd2=zeros(n2,r2);
            uadd3=zeros(n3,r3);
            %Compute ``maximum volume'' indices from sbm
            ind1=maxvol2(u1,ind1);
            ind2=maxvol2(u2,ind2);
            ind3=maxvol2(u3,ind3);
            [ind01,ind02,ind03]=maxvol3d(core);
            
            %Extract true rows
            if ( iter < niters )
                for i1=1:n1
                    for k=1:r1
                        uadd1(i1,k)=f(i1,ind1(ind01(k,1)),ind2(ind01(k,2)),mat);
                        
                    end
                end
                for i2=1:n2
                    for k=1:r2
                        uadd2(i2,k)=f(ind1(ind02(k,2)),i2,ind3(ind02(k,1)),mat);
                    end
                end
                for i3=1:n3
                    for k=1:r3
                        uadd3(i3,k)=f(ind1(ind03(k,1)),ind2(ind03(k,2)),i3,mat);
                    end
                end
            end
        end
        return
    end
    function [b]=convt(a,u1,u2,u3)
        [n1,n2,n3]=size(a);
        [n1,r1]=size(u1);
        [n2,r2]=size(u2);
        [n3,r3]=size(u3);
        b=a;
        b=u3'*(reshape(b,[n1*n2 n3]))';
        b=u2'*(reshape(b,[r3*n1 n2]))';
        b=u1'*(reshape(b,[r2*r3 n1]))';
        b=reshape(b,[r1,r2,r3]);
        return
    end

% a is n x r matrix, we have to find r rows that give a good submatrix
    function [ind]=maxvol2(a,ind0)
        [n,r]=size(a);
        if ( n == r )
            ind = 1:r;
            return
        end
        %ind=zeros(r,1);
        %Initialize
        [l,u,p]=lu(a,'vector');
        ind=p(1:r);
        %ind=1:r;
        %p=1:n;
        %ind=ind0;
        %p=1:n;
        %for i=1:r
        %    p(i)=ind(i);
        %    p(ind(i))=i;
        %end
        b=a(p,:);
        sbm=a(ind,:);
        z=b(r+1:n,:)*inv(sbm);
        %Start iterations
        niters =20;
        gm=2;
        eps=1e-4;
        iter=0;
        while (iter <= niters);
            [mx,pv]=max((abs(z))');
            [mx0,pv0]=max(mx);
            if ( mx0 <= 1 + eps )
                break;
            end
            %Swap i0-th row with j0-th row
            i0=pv0;
            j0=pv(i0);
            i=p(i0+r);
            p(i0+r)=p(j0);
            p(j0)=i;
            bb=z(:,j0);
            bb(i0)=bb(i0)+1;
            cc=z(i0,:);
            cc(j0)=cc(j0)-1;
            z=z-bb*cc./z(i0,j0);
            iter=iter+1;
            ind=p(1:r);
            %fprintf('%e, %e \n',abs(det(b(ind,:))),abs(det(b(ind,:))*mx0))
            %fprintf('log(det): %f\n',log(abs(det(b(ind,:)))));
            %please check consistency
        end
        %Dummy check;
        %a0=abs(a*inv(a(ind,:)));
        %fprintf('Maximal element: %f \n',max(abs(a0(:))));
        return
    end

    function [ind1,ind2,ind3]=maxvol3d(a)
        [n1,n2,n3]=size(a);
        %Compute unfoldings
        a3=reshape(a,[n1*n2,n3]);
        a2=reshape(a3',[n3*n1,n2]);
        a1=reshape(a2',[n2*n3,n1]);
        ind1=1:n1;
        ind2=1:n2;
        ind3=1:n3;
        ind1=maxvol2(a1,ind1);
        ind2=maxvol2(a2,ind2);
        ind3=maxvol2(a3,ind3);
        ind01=zeros(n1,2);
        ind02=zeros(n2,2);
        ind03=zeros(n3,2);
        for i1=1:n1
            ind=ind1(i1);
            i3=floor((ind-1)/n2)+1;
            i2=ind-(i3-1)*n2;
            ind01(i1,1)=i2;
            ind01(i1,2)=i3;
        end
        for i2=1:n2
            ind=ind2(i2);
            i1=floor((ind-1)/n3)+1;
            i3=ind-(i1-1)*n3;
            ind02(i2,1)=i3;
            ind02(i2,2)=i1;
        end
        for i3=1:n3
            ind=ind3(i3);
            i2=floor((ind-1)/n1)+1;
            i1=ind-(i2-1)*n1;
            ind03(i3,1)=i1;
            ind03(i3,2)=i2;
        end
        ind1=ind01;
        ind2=ind02;
        ind3=ind03;
        return
    end

    function [r]=my_chop(a,eps)
        nrm2=norm(a).^2;
        bound=eps.^2*nrm2;
        [r1,tmp]=size(a(:));
        er=0;
        r=r1;
        while(er < bound && r >= 1 )
            er=er+a(r).^2;
            r=r-1;
        end
        r=r+1;
        return
    end


    function [u1,u2,u3,core]=svd3l(a,eps)
        [n1,n2,n3]=size(a);
        %Compute unfoldings
        a3=reshape(a,[n1*n2,n3]);
        a2=reshape(a3',[n3*n1,n2]);
        a1=reshape(a2',[n2*n3,n1]);
        [tmp,s3,u3]=svd(a3,'econ');
        s3=diag(s3);
        r3=my_chop(s3,eps/3);
        [tmp,s2,u2]=svd(a2,'econ');
        s2=diag(s2);
        r2=my_chop(s2,eps/3);
        [tmp,s1,u1]=svd(a1,'econ');
        s1=diag(s1);
        r1=my_chop(s1,eps/3);
        r1=min(r1,n1);
        r2=min(r2,n2);
        r3=min(r3,n3);
        u1=u1(:,1:r1);
        u2=u2(:,1:r2);
        u3=u3(:,1:r3);
        core=convt(a,u1,u2,u3);
        return
    end


%% end
end