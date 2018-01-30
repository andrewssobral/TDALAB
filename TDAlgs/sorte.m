function [n ssig e]=sorte(x)
%% Implementation of Second ORder sTatistic of the Eigenvalues
% (SORTE)
[M T]=size(x);
if M>T
    x=x';
    [M T]=size(x);
end
if M>1
    [u s]=eig(cov(x'));
    [s pos]=sort(diag(s),'descend');
    T=length(s);
else
    s=sort(x,'descend');
end
e=s;

if T<=5
    n=sum(s>s(1)*0.05);
    ssig=e;
    return;
end

ds=-(diff(s));
for k=1:T-1
    sig(k)=mean((ds(k:T-1)-mean(ds(k:T-1))).^2);
    if k>1
        ssig(k-1)=sig(k)./max(sig(k-1));
    end
end
[K n]=min(ssig(1:T-4));
ssig(T-3:end)=[];
