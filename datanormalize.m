function [x xnorm]=datanormalize(x,nor,ori)
% Usage: function x=normalize(x,nor,ori);
%      x=normalize(x) for short.
% normalize the matrix x by specified norm(nor) and orientation (ori)
% nor: {2|1|p|inf} norm used. [default 2]
% ori: only valid for a matrix input.
%      1,'col' -- column-wise normalization [default]
%      2,'row' -- row-wise normalization
% Coded by Guoxu Zhou
if nargin==1
    nor=[];    ori=[];
elseif nargin==2
    ori=[];
end
if isempty(nor) nor=2; end
if isempty(ori) ori='col'; end
switch ori
    case {'col',1}
        ori=1;
    otherwise
        ori=2;
end

typeofx=class(x);
switch typeofx
    case 'double'
        if numel(size(x))~=2
            error('For double data, a matrix is expected.');
        end
        nor=max(nor,0.1);
        if nor==2
            xnorm=max(sum(x.*x,ori),eps).^0.5;
        elseif nor==1
            xnorm=max(sum(abs(x),ori),eps);
        elseif isinf(nor)
            xnorm=max(max(abs(x),[],ori),eps);
        else
            xnorm=max(sum(abs(x).^nor,ori),eps).^(1/nor);
        end
        x=bsxfun(@rdivide,x,xnorm);
    case 'ktensor'
        N=numel(size(x));
        for n=1:N
            [x.U{n} nu]=datanormalize(x.U{n},nor,'col');
            x.lambda=x.lambda.*nu';
        end
    case 'ttensor'
        N=numel(size(x));
        for n=1:N
            [x.U{n} nu]=datanormalize(x.U{n},nor,'col');
            x.core=ttm(x.core,diag(nu),n);
        end
    otherwise
        error('Unsuported data type.');
end
end
        
        
