% #########################################################################
%  Copyright C. Caiafa and A. Cichocki 2010.
%  Laboratory for Advanced Brain Signal Processing
% Written by C. Caiafa, 2010.
% 
% Fiber Sampling Tensor Decomposition type 2 (FSTD2) for 3D-tensors(N=3) as
% an application of the FSTD2 formula presented in "Generalizing the Column-Row Matrix Decomposition to Multi-way Arrays" by
% C. Caiafa and C. Cichocki, Linear Algebra and its Applications, Vol. 433,
% pp. 557?73, 2010 (Elsevier). The fibers and indices are NOT selected in an optimal
% way (are just by random). We compare FSTD2 to the approximation obtained by the classical Tucker-ALS algorithm. 
% This program requires to have the Tensor Toolbox installed: Brett W. Bader and Tamara G. Kolda, MATLAB Tensor Toolbox Version 2.4, http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox/, March 2010.

function [Yaprox FIB U fit]=FSTDT2rnd(Y,opts)
% Y: tensor to be decomposed
% k: Number of indices to select in each mode
defopts = struct('Rank',[],'Verbose',false);
if ~exist('opts','var')
    opts = struct;
end
opts=scanparam(defopts,opts);

%% Set algorithm parameters from input or by using defaults
k=opts.Rank;
Y=double(Y);
normY=Y(:);normY=sqrt(normY'*normY);
I=size(Y); % mode sizes
if length(k)<length(I)
    k(end+1:length(I))=k(end);
end

N=size(I,2); % dimensions

indexfib=cell(1,N); % indices of fibers in the corresponding unfolding matrix
index=cell(1,N); % indices of subtensor
FIB=cell(1,N); %selected fibers
Wsub=cell(1,N); %selected sub matrices within fibers
p=prod(I);
for n=1:N
    % random selection of fibers
    randindex=randperm(p/I(n));    
    indexfib{n}=randindex(1:k(n));
    Yaux=double(tenmat(tensor(Y),n));
    FIB{n}=Yaux(:,indexfib{n}); 
    % random selection of subtensor
    randindex=randperm(I(n));    
    index{n}=randindex(1:k(n));
    Wsub{n}=FIB{n}(index{n},:);
    Wsub{n}=inv(Wsub{n});
end
W=Y(index{:});
 
U=ttensor(tensor(W),Wsub);
Yaprox=ttensor(U,FIB{:}); 


z=double(Yaprox)-Y;z=z(:);
fit=1-sqrt(z'*z)/normY;

sum=0;
p=prod(k);
for n=1:N
    sum=sum+I(n)*p/k(n);
end
ratio=(sum-(N-1)*p)/prod(I);
if opts.Verbose
    fprintf('Fitting: %f   Sampling ratio: %f\n',fit,   ratio);
end
end
