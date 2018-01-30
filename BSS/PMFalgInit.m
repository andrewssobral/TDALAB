function [PMFparaList PMFalgs PMFparaTypeList]=PMFalgInit()
%% penalty matrix factorization algs definition.

global algs paraList paraTypeList;

allidx=1:numel(algs);
flag=[algs(:).model]==AlgType.PMF;
allidx=allidx(flag);

PMFparaList=paraList(allidx);
PMFalgs=algs(allidx);
PMFparaTypeList=paraTypeList(allidx);





















%% Add your algorithm BEFORE this line.
%% ----------------------------------------------------
%% END OF THIS FILE. DON'T CHANGE ANYTHING AFTER HERE.

%%%%%%%%%%%%%%%%%%%%%%%
% ## paraType:
% paraType.paraTF
% paraType.paraMenu
% paraType.paraDouble
% paraType.paraString

%%
if nargout>1
    allidx=1:numel(algs);
    
end