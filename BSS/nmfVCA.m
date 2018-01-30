function [w h]=nmfVCA(x,opts)
% Vertex Component Analysis Algorithm [VCA]
%
% [VCA] J. Nascimento and J. Bioucas-Dias
% "Vertex component analysis: a fast algorithm to unmix hyperspectral data"
% IEEE Transactions on Geoscience and Remote Sensing,  vol. 43, no. 4, 
% pp. 898-910, 2005.
%
% -------------------------------------------------------------------
% Usage:
%
% [Ae, indice, Rp ]= vca(R,'Endmembers',p,'SNR',r,'verbose',v)
%
% ------- Input variables -------------------------------------------
%
%  R - matrix with dimensions L(channels) x N(pixels)
%      Each pixel is a linear mixture of p endmembers
%      signatures R = M X, where M  and X are the mixing matrix 
%      and the abundance fractions matrix, respectively.
%
% 'Endmembers'
%          p - number of endmembers in the scene
%
% ------- Output variables -------------------------------------------
%
% A      - estimated mixing matrix (endmembers signatures)
%
% indice - pixels chosen to be the most pure
%
% Rp     - Data R projected on the identified signal subspace
%
% ------- Optional parameters -----------------------------------------
%
% 'SNR'     - (double) signal to noise ratio (dB)
%             SNR is used to decide the type of projection: projective
%             or orthogonal.
%
% 'verbose' - [{'on'} | 'off']
% ---------------------------------------------------------------------
%
% Please see [VCA] for more details or contact the Authors
%
% -----------------------------------------------------------------------
% version: 3.0 (21-January-2012)
%
% Modifications w.r.t. version 2.1:
%     
%  - Increased efficiency in the memory usage
%  - Correction of a bug in SNR estimation
%  - detection of outliers in the projective projection
%
% -----------------------------------------------------------------------
% Copyright (2012):  Jos?Nascimento (zen@isel.pt)
%                    Jos?Bioucas Dias (bioucas@lx.it.pt)
%
% affineProj is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."

%% FOR TDALAB 
% X is transposed
x=x';

defopts=struct('NumOfComp',[],'verbose',false,'SNR',[]);
if ~exist('opts','var')
    opts=struct();
end
[r,snr,verbose]=scanparam(defopts,opts);
[M T]=size(x);
if isempty(r)
    r=M;
end

x=max(x,0);
x0=x;
sn=sum(x);
flag=sn<eps;
sn(flag)=[];
x(:,flag)=[];

x=bsxfun(@rdivide,x,sn);
if isempty(r)
    w=VCA(x,'Endmembers',r,'verbose',verbose);
else
    w=VCA(x,'Endmembers',r,'SNR',snr,'verbose',verbose);
end
h=w';
w=max(w\x0,0)';