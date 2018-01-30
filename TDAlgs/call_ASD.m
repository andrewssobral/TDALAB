function [P res]=call_ASD(Y,opts)
% Usage: [P res]=callASD(Y,opts)
%  Opts: [NumOfComp],[Tol],[lambda],[MaxIter]
%
% Core codes are provided by the original author and this
% version is enclosed for calling in TDALAB.
% Ref:
% Three-way data resolution by alternating slice-wise diagonalization (ASD) method
defopts = struct('NumOfComp',0,'Tol',1e-10,'lambda',1e-3,'MaxIter',10000);
if ~exist('opts','var')
    opts = struct;
end
[J,tol,lambda,maxIter]=scanparam(defopts,opts);

if strcmp(class(Y),'ktensor')||strcmp(class(Y),'ttensor')||strcmp(class(Y),'tensor')
    Y=double(Y);
end


if iscell(Y)
    if length(Y)~=3
        error('only works for PARAFAC/CP-3D model.');
    end
    
    [II NN]=size(Y{1});
    JJ=size(Y{2},1);
    KK=size(Y{3},1);
    rrzr=[];rrzl=[];
    for k=1:KK
        rrzr=[rrzr Y{1}*diag(Y{3}(k,:))*Y{2}'];
        rrzl=[rrzl Y{2}*diag(Y{3}(k,:))*Y{1}'];
    end
elseif strcmp(class(Y),'double')&&ndims(Y)==3
    KK=size(Y,3);
    rrzr=[];rrzl=[];
    for k=1:KK
        rrzr=[rrzr Y(:,:,k)];
        rrzl=[rrzl Y(:,:,k)'];
    end    
else
    error('Unsuported input data type.');
end


[Z, res] = ASD(rrzr, rrzl, J,tol,lambda,maxIter);
P=ktensor(Z(:));
%%=========================================================================
% Notations:
% rrzk is the response matrix of the kth sample.
% KK is the number of samples. KK > = 2. In the study we use KK = 4. --
% length of z
% rrzr = [rrz1 rrz2 ?rrzKK];
% rrzl = [rrz1' rrz2' ?rrzKK'];
% II is the number of variables in x order. In the study II = 50.
% JJ is the number of variables in y order. In the study JJ = 20.
% XX, YY and ZZ are the resolved profiles in x, y and z orders respectively.
% NN is an estimate of the number of components.  -- rank
% uux and uuy are the first N singular vectors of rrzr * rrzr' and rrzl * rrzl' respectively.
% weight is the penalty weight.  -- \lambda
% AA = uux'* XX; BB = uuy'* YY; GG = inv(AA'); HH = inv(BB').
function [A, errflow] = ASD(rrzr, rrzl, NN,tol,weight,maxIter)
II=size(rrzr,1);JJ=size(rrzl,1);
KK=size(rrzr,2)/JJ;


% II is unused.
[rrzr2, rrzl2, uux, uuy] = svddr(rrzr, rrzl, NN, JJ, KK);
GG = eye(NN); HH = eye(NN);
err1 = 100; derr = 1; cyc = 0;;
errflow = [];
while derr > tol & cyc < maxIter
    weight=weight*1.2;
    cyc = cyc + 1;
    err0 = err1;
    ZZ = gghhtozz(rrzr2, GG, HH, NN, KK);
    AA = inv(GG'); BB = inv(HH');
    err1 = error1(rrzr2, GG, HH, KK, NN);
    HH = ggzztohh(rrzr2, GG, ZZ, BB, weight, NN, KK);
    HH = HH * diag(ones(1, NN) ./ sqrt(sum(HH .* HH)));
    ZZ = gghhtozz(rrzr2, GG, HH, NN, KK);
    GG = ggzztohh(rrzl2, HH, ZZ, AA, weight, NN, KK);
    GG = GG * diag(ones(1, NN) ./ sqrt(sum(GG.* GG)));
    errflow = [errflow, err1];
    derr = abs(err0 - err1);
end
XX = uux * AA; YY = uuy * BB;
XX = XX * diag(ones(1, NN) ./ sqrt(sum(XX.* XX)));
YY = YY * diag(ones(1, NN) ./ sqrt(sum(YY.* YY)));
ZZ = xxyytozz(rrzr, XX, YY, JJ, KK, NN);
% Post-processing to keep sign convention.
[maxx, indx] = max(abs(XX));
[maxy, indy] = max(abs(YY));
xsign = ones(NN, 1); ysign = ones(NN, 1);
for nn = 1:NN
    xsign(nn) = sign(XX(indx(nn), nn));
    ysign(nn) = sign(YY(indy(nn), nn));
end
XX = XX * diag(xsign);
YY = YY * diag(ysign);
ZZ = ZZ * diag(xsign) * diag(ysign);

A{1}=XX;A{2}=YY;A{3}=ZZ;

function err1 = error1(rrzr, pp, qq, KK, NN)
drrz = [];
for kk = 1:KK
    rrz = pp' * rrzr(:, NN * (kk - 1) + 1:NN * kk) * qq;
    drrz = [drrz, rrz - diag(diag(rrz))];
end
err1 = sum(sum(drrz.*drrz));

function ZZ = gghhtozz (rrzr, GG, HH, NN, KK)
for kk = 1:KK
    ZZ(kk, :) = diag (GG' * rrzr (:, NN * (kk-1) + 1:NN * kk) * HH)';
end

function HH = ggzztohh (rrzr, GG, ZZ, BB, weight, NN, KK)
tt1 = zeros (NN, NN);
tt2 = zeros (NN, NN);
for kk = 1:KK
    tt1 = tt1 + rrzr (:, NN * (kk-1) + 1:NN * kk)' * GG * GG' * rrzr (:, NN * (kk-1) + 1:NN * kk);
    tt2 = tt2 + rrzr (:, NN * (kk-1) + 1:NN * kk)' * GG * diag(ZZ(kk,:));
end
HH = inv (tt1 + weight * BB * BB') * (tt2 + weight * BB);

function [rrzr2, rrzl2, uux, uuy] = svddr (rrzr, rrzl, NN, JJ, KK)
[uu, ss, vv] = svd (rrzl * rrzl');
uuy = uu (:, 1:NN);
[uu, ss, vv] = svd (rrzr * rrzr');
uux = uu (:, 1:NN);
rrzr2 = []; rrzl2 = [];

%% HERE
for kk = 1:KK
    rrzr2 = [rrzr2 uux' * rrzr(:, JJ * (kk-1) + 1:JJ * kk) * uuy];
    rrzl2 = [rrzl2 uuy' * rrzr(:, JJ * (kk-1) + 1:JJ * kk)' * uux];
end

function ZZ = xxyytozz (rrzr, XX, YY, JJ, KK, NN)
FF = zeros (KK, NN);
for kk = 1:KK
    FF(kk,:) = diag (XX' * rrzr (:, JJ * (kk-1) + 1:JJ * kk) * YY)';
end
ZZ = FF * inv ((XX' * XX).* (YY' * YY));

