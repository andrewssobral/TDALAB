function [SIR maps] = CalcSIR(A,Aest)
% Sergio Cruces & Andrzej Cichocki
% A(:,maps) matches Aest properly.  --Added by Guoxu Zhou

% mean value should be extracted first --Added by Guoxu Zhou
A=real(A);
A=bsxfun(@minus,A,mean(A));
Aest=bsxfun(@minus,Aest,mean(Aest));

A=A*diag(1./(sqrt(sum(A.^2))+eps));
Aest=Aest*diag(1./(sqrt(sum(Aest.^2))+eps));

col=size(A,2);flag=zeros(1,col);
MSE=inf*ones(1,col);
for i=1:size(Aest,2)
    temp=min(sum(bsxfun(@minus,Aest(:,i),A).^2,1),...
        sum(bsxfun(@plus,Aest(:,i),A).^2,1));
    temp=max(temp,flag);
    [MSE(i),maps(i)]=min(temp);
    flag(maps(i))=inf;
end
SIR=-10*log10(MSE);


