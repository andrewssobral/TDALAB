function Y=NNTuckerCore(Y,U)
%% Given tensor/ttensor Y
% find the optimal nonnegative core tensor with fixed nonnegative loading matrices U.
% If U is not given, U=Y.U.

    if nargin==1
        if strcmp(class(Y),'ttensor')
            U=Y.U;
            R=size(Y.core);
        else
            error('U is required.');
        end
    end
    nc=rand(R);
    enum=max(double(ttm(Y,U,'t')),eps);
    UUT=cellfun(@(x) x'*x,U,'uni',false);
    for i=1:200
        nc=nc.*(enum./max(double(ttm(tensor(nc),UUT)),eps));
    end
    
    Y=ttensor(tensor(nc),U);
end

