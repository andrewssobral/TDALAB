function [ Ycap ] = n3cp( Y,opts )
%% N-way CP decomposition based on a proceeding 3-way decomposition
%  Usage: Ycap=n3cp(Y,opts);
%  opts
%      .NumOfComp: rank of Ycap
%      .modeGroup: group of the modes. Such as: 1 2 ; 3 ; 4 of a 4-way
%          tensor. Split by ';'.
%      .CP3_FILE: algorithm for 3-way CP
%      .lra: none|pca|sampling low-rank approximation method
%      .lra_rank: rank of lra (i.e., reduced dimension)
%      .krp_init: random|svd initialization method for Khatri-rao prod
%           approximation
%      .krp_maxit: maxiteration 
%      .krp_costol: stopping criterion 
%      .krp_constraints: none|nonnegative|sparse
%      .krp_cons_para: parameter which is used for krp_constraints.
