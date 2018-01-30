function [W H]=PMFSmoISNMF(V,opts)
%% Smooth Itakura-Saito NMF solving
%
% min D(V|WH) + lambda * sum_kn d( h_k(n-1) | h_kn )
%
% [W, H, cost] = nmf_is_smooth(V, n_iter, W_ini, H_ini, lambda)
%
% Input:
%   - V: positive matrix data (F x N)
%   - n_iter: number of iterations
%   - W_ini: dictionary init. (F x K)
%   - H_ini: activation coefficients init. (K x N)
%   - lambda: penalty weight (positive scalar)
% Outputs :
%   - W and H such that
%
%               V \approx W * H
%
%   - cost: IS divergence though iterations
%   - cost_pen: penalized cost function
%
% Reference paper
%
% C. Fevotte, "Majorization-minimization algorithm for smooth Itakura-Saito
% nonnegative matrix factorizations", Proc. ICASSP 2011.
%
% Copyright (C) 2010 C. Fevotte
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  @@ Interface Modified by Guoxu Zhou for the purpose of TDALAB @@
%
defopts=struct('NumOfComp',[],'n_iter',200,'lambda',10);
if ~exist('opts','var')
    opts=struct;
end
[NComps,n_iter,lambda]=scanparam(defopts,opts);

% function [W, H, cost, cost_pen] = nmf_is_smooth(V, n_iter, W, H, lambda)

%
% Outputs :
%   - W and H such that
%
%               V \approx W * H
%
%   - cost: IS divergence though iterations
%   - cost_pen: penalized cost function
%
% Reference paper
%
% C. Fevotte, "Majorization-minimization algorithm for smooth Itakura-Saito
% nonnegative matrix factorizations", Proc. ICASSP 2011.
%
% Copyright (C) 2010 C. Fevotte
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



[F,N] = size(V);

%% Initialization for NC, W, H.
if isempty(NComps)
    NComps=min(F,N);
end
W=rand(F,NComps);
H=rand(NComps,N);
% !! End of initialization

K = size(W,2);

cost = zeros(1,n_iter);
cost_pen = zeros(1,n_iter);

% Compute data approximate
V_ap = W*H;

% Compute initial cost value
cost(1) = sum(V(:)./V_ap(:) - log(V(:)./V_ap(:))) - F*N;

R = H(:,1:N-1)./H(:,2:N);
cost_pen(1) = cost(1) + lambda * (sum(R(:)-log(R(:)))-K*(N-1));

for iter=2:n_iter
    
    %% Update H %%
    Gp = W'*V_ap.^-1;
    Gn = W'*(V.*V_ap.^-2);
    Ht = H;
    
    % n = 1 %
    p2 = Gp(:,1) + lambda./H(:,2);
    p1 = - lambda;
    p0 = - Gn(:,1).*Ht(:,1).^2;
    H(:,1) = (sqrt(p1.^2-4*p2.*p0) - p1)./(2*p2);
    
    % n = 2:N-1 %
    for n=2:N-1
        H(:,n) = sqrt((Gn(:,n).*Ht(:,n).^2 + lambda*H(:,n-1))./(Gp(:,n) + lambda./H(:,n+1)));
    end
    
    % n = N %
    p2 = Gp(:,N);
    p1 = lambda;
    p0 = - (Gn(:,N).*Ht(:,N).^2 + lambda*H(:,N-1));
    H(:,N) = (sqrt(p1.^2-4*p2.*p0) - p1)./(2*p2);
    
    V_ap = W*H;
    
    %% Update W %%
    W = W .* ((V.*V_ap.^-2)*H')./(V_ap.^-1*H');
    V_ap = W*H;
    
    %% Normalization %%
    scale = sum(W,1);
    W = W .* repmat(scale.^-1,F,1);
    H = H .* repmat(scale',1,N);
    
    % Compute costs
    cost(iter) = sum(sum(V./V_ap - log(V./V_ap))) - F*N;
    
    R = H(:,1:N-1)./H(:,2:N);
    cost_pen(iter) = cost(iter) + lambda * (sum(R(:)-log(R(:)))-K*(N-1));
    
end

