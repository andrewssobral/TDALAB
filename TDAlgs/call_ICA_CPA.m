function [ Ycap] = call_ICA_CPA( Y,opts)
%% This function calls the ICA_CPA algorithm.
% Usage: Ycap=call_ICA_CPA(Y,opts);
% Output: Ycap is a ktensor (CP tensor)
% Input: Y -- tensor
%        opts.
%            R: rank of Ycap;
%            alfa: weighting factor (between 0 and pi/2) that weighs the relative contribution of
%                   second order and higher order statistics, low value implies high
%                   contribution of higher order statistics
%
% The original information about this algorithm is below:
%function [S2,A2,B2] = ica_cpa (tensor,R,alfa)
%
%This function estimates the trilinear CP decomposition from a tensor in
%which one mode contains statistically independent sources.
%
% INPUT: tensor: 3 dimensional tensor with independence assumed in third
% dimension
%        R : number of sources to be estimated                                                                                                
%        alfa: weighting factor (between 0 and pi/2) that weighs the relative contribution of
%        second order and higher order statistics, low value implies high
%        contribution of higher order statistics
%
% OUTPUT: mode 3: S2 : statistically independent sources
%         mode 1: A2
%         mode 2: B2
%
%
%De Vos M., Nion D., Van Huffel S. and De Lathauwer L., 
%"A combination of Parallel Factor and Independent Component Analysis", 
%Signal Processing, vol. 92, Jul. 2012, pp. 2990-2999.
%--------------------------------------------------------------------------------------------------------
% @Copyright 2010    
% Lieven De Lathauwer, Maarten De Vos, Dimitri Nion 
% lieven.delathauwer@kuleuven-kortrijk.be
% K.U. Leuven, campus Kortrijk, Belgium
% This M-file and the code in it belongs to the holder of the copyrights. 
%-------------------------------------------------------------------------------------------------------

defopts=struct('NumOfComp',[],'alfa',[]);
if ~exist('opts','var')
    opts=struct();
end
[R alfa]=scanparam(defopts,opts);
Y=double(Y);
[S2,A2,B2] = ica_cpa(Y,R,alfa);
Ycap=ktensor({A2,B2,S2});

end

