function [W, H, cost] = beta_nmf_ME(V,opts)
% Majorization-Equalizaton algorithm for NMF with the beta-divergence
%
% Implemented for beta = {-1,0,0.5,1.5,2,3}
%
% [W, H, cost] = beta_nmf_ME(V, beta, n_iter, W_ini, H_ini)
%
% Input:
%   - V: positive matrix data (F x N)
%   - beta : beta-divergence parameter {-1,0,0.5,1.5,2,3}
%   - n_iter: number of iterations
%   - W_ini: basis (F x K)
%   - H_ini: gains (K x N)
%
% Outputs :
%   - W and H such that
%
%               V \approx W * H
%
%   - cost : beta-divergence though iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2011 Cedric Fevotte & Jerome Idier
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% Please cite the following reference if you use this code
% C. Fevotte & J. Idier, "Algorithms for nonnegative matrix factorization
% with the beta-divergence ", Neural Compuation, 2011.
%
% Please report any bug at fevotte -at- telecom-paristech -dot- fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


defopts=struct('NumOfComp',[],'beta',0.5,'n_iter',200);
if ~exist('opts','var')
    opts=struct;
end
[NComps,beta,n_iter]=scanparam(defopts,opts);

% function [W, H, cost] = beta_nmf_ME(V, beta, n_iter, W, H)

% Majorization-Equalizaton algorithm for NMF with the beta-divergence
%
% Implemented for beta = {-1,0,0.5,1.5,2,3}
%
% [W, H, cost] = beta_nmf_ME(V, beta, n_iter, W_ini, H_ini)
%
% Input:
%   - V: positive matrix data (F x N)
%   - beta : beta-divergence parameter {-1,0,0.5,1.5,2,3}
%   - n_iter: number of iterations
%   - W_ini: basis (F x K)
%   - H_ini: gains (K x N)
%
% Outputs :
%   - W and H such that
%
%               V \approx W * H
%
%   - cost : beta-divergence though iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2011 Cedric Fevotte & Jerome Idier
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% Please cite the following reference if you use this code
% C. Fevotte & J. Idier, "Algorithms for nonnegative matrix factorization
% with the beta-divergence ", Neural Compuation, 2011.
%
% Please report any bug at fevotte -at- telecom-paristech -dot- fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch_W = 1; 
switch_H = 1;
SCALE = and(switch_W,switch_H) & 1;

theta = 0.95; % Relaxation coefficient

[F,N] = size(V);

%% initialization for NComps, W, H
if isempty(NComps)
    NComps=min(F,N);
end
W=rand(F,NComps);
H=rand(NComps,N);
% !! End of initialization

V_ap = W*H;

cost = zeros(1,n_iter);
cost(1) = betadiv(V,V_ap,beta);


switch beta
    
    case -1 % ME update always defined
        
        for iter = 2:n_iter
            
            if switch_W
                W_h = W .* ((V.*V_ap.^-3)*H')./(V_ap.^-2*H');
                W = ((W_h.^2+8*W_h.*W).^.5+W_h)/4;
                V_ap = W*H;
            end
            
            if switch_H
                H_h = H .* (W'*(V.*V_ap.^-3))./(W'*V_ap.^-2);
                H = ((H_h.^2+8*H_h.*H).^.5+H_h)/4;
                V_ap = W*H;
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,V_ap,beta);
            
        end
        
    case 0 % ME update always defined
                
        for iter = 2:n_iter
            
            if switch_W
                W = W .* ((V.*V_ap.^-2)*H')./(V_ap.^-1*H');
                V_ap = W*H;
            end
            
            if switch_H
                H = H .* (W'*(V.*V_ap.^-2))./(W'*V_ap.^-1);
                V_ap = W*H;
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,V_ap,beta);
            
        end
        
    case .5 % ME update always defined

        for iter = 2:n_iter
            
            if switch_W
                W_h = W .* ((V.*V_ap.^-1.5)*H')./(V_ap.^-0.5*H');
                W = (sqrt(W+8*W_h)-sqrt(W)).^2/4;
                V_ap = W*H;
            end
            
            if switch_H
                H_h = H .* (W'*(V.*V_ap.^-1.5))./(W'*V_ap.^-0.5);
                H = (sqrt(H+8*H_h)-sqrt(H)).^2/4;
                V_ap = W*H;
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,V_ap,beta);
            
        end
        
    case 1.5 % mixed update prolonged ME and MM
        
        for iter = 2:n_iter
                        
            if switch_W                
                W_h = W .* ((V.*V_ap.^-0.5)*H')./(V_ap.^0.5*H');                
                DW = 12*W_h-3*W;
                IW = 3*W_h<W;
                W = real((sqrt(DW)-sqrt(W)).^2/4);
                W(IW) = 0;
                W = W_h+theta*(W-W_h);
                V_ap = W*H;                
            end
            
            if switch_H                
                H_h = H .* (W'*(V.*V_ap.^-0.5))./(W'*V_ap.^0.5);                
                DH = 12*H_h-3*H;
                IH = 3*H_h<H;
                H = real((sqrt(DH)-sqrt(H)).^2/4);
                H(IH) = 0;
                H = H_h+theta*(H-H_h);
                V_ap = W*H;
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,V_ap,beta);
            
        end
                
    case 2 % mixed update of prolonged ME and MM
        
        for iter = 2:n_iter
            
            if switch_W                
                W_h = W .* (V*H')./(W*(H*H'));
                W = 2*W_h-W;
                W(W<=0) = 0;
                W = W_h+theta*(W-W_h);
            end
            
            if switch_H                
                H_h = H .* (W'*V)./((W'*W)*H);
                H = 2*H_h-H;
                H(H<=0) = 0;
                H = H_h+theta*(H-H_h);
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,W*H,beta);
            
        end
        
    case 3 % mixed update of prolonged ME and heuristic
        
        for iter = 2:n_iter
            
            if switch_W               
                W_h = W .* ((V.*V_ap)*H')./(V_ap.^2*H');
                DW = W.*(12*W_h-3*W);
                IW = 3*W_h<W; 
                W = real(sqrt(DW)-W)/2;
                W(IW) = 0;
                W = W_h+theta*(W-W_h);
                V_ap = W*H;
            end
            
            if switch_H                
                H_h = H .* (W'*(V.*V_ap))./(W'*V_ap.^2);
                DH = H.*(4*beta*H_h-3*H);
                IH = 3*H_h<H;
                H = real(sqrt(DH)-H)/2;
                H(IH) = 0;
                H = H_h+theta*(H-H_h);
                V_ap = W*H;
            end
            
            if SCALE
                scale = sum(W,1);
                W = W .* repmat(scale.^-1,F,1);
                H = H .* repmat(scale',1,N);
            end
            
            cost(iter) = betadiv(V,V_ap,beta);
            
        end
        
end


%% AUX func
function d = betadiv(A,B,beta)

% beta-divergence
%
%   d = betadiv(A,B,beta)
% 
% - a \= 0,1
%   d(x|y) = ( x^a + (a-1)*y^a - a*x*y^(a-1) ) / (a*(a-1))
% - a = 1
%   d(x|y) = x(log(x)-log(y)) + (y-x)  KULLBACK-LEIBLER
% - a = 0
%   d(x|y) = x/y - log(x/y) -1         ITAKURA-SAITO

switch beta
    case 2
        d = sum((A(:)-B(:)).^2)/2;
    case 1
        ind_0 = find(A(:)<=eps);
        ind_1 = 1:length(A(:));
        ind_1(ind_0) = [];
        d = sum( A(ind_1).*log(A(ind_1)./B(ind_1)) - A(ind_1) + B(ind_1) ) + sum(B(ind_0));        
    case 0
        d = sum( A(:)./B(:) - log(A(:)./B(:)) ) - length(A(:));
    otherwise
        d = sum( A(:).^beta + (beta-1)*B(:).^beta - beta*A(:).*B(:).^(beta-1) )/(beta*(beta-1));
end


