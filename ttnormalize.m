function [ Y ] = ttnormalize( Y,nrm )
if nargin==1
    nrm=[];
end
if isempty(nrm)
    nrm=0;
end

if strcmpi(class(Y),'ttensor')
    I=size(Y);
    N=numel(I);
    J=size(Y.core);
    switch nrm
        case 0
            for n=1:N
                ts=0:I(n):(J(n)-1)*I(n);
                [temp idx]=max(abs(Y.U{n}));
                d=Y.U{n}(idx+ts);
                d=sign(d).*max(abs(d),eps);
                Y.U{n}=bsxfun(@rdivide,Y.U{n},d);
                Y.core=ttm(Y.core,diag(d),n);
            end
        case 1
            for n=1:N
                d=max(sum(abs(Y.U{n})),eps);
                Y.U{n}=bsxfun(@rdivide,Y.U{n},d);
                Y.core=ttm(Y.core,diag(d),n);
            end                
        case 2
            for n=1:N
                d=max(sum(Y.U{n}.^2),eps);
                Y.U{n}=bsxfun(@rdivide,Y.U{n},d);
                Y.core=ttm(Y.core,diag(d),n);
            end
        otherwise
            disp('Unsupported norm.');
    end
    

end

end

