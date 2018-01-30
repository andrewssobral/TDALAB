function [str]=defstr(type)
type=lower(type);
switch type
    %% For input format
    case 'input'
        str={'ktensor','ttensor'};
        
        %% Tensor decomposition model str
    case 'tdmodel'
%         str={'CP','Tucker','Parafac2','BCD','ParaTucker','DEDICOM'};
        str={'CP','Tucker','PMF','BCD'};
    
        %% Tucker mode constraint type string
    case 'tuckerconstr'
        str={'None','Penalized Matrix Factorization'}; %% Must no more than five
        
        %% Algorithms
    case 'penalized matrix factorization'
        [paraList PMFalgs]=PMFalgInit;
        str={'Free',PMFalgs(:).details};
        
        %% core tensor
    case 'gencore' 
        str={'Auto';'Uniform';'Gaussian';'Abs. Gaussian'};
        
        %% ortho
%     case 'orthogonality'
%         str={'Free','orth'};
%         
%         %% SmoCA
%     case 'smoothness'
%         str={'Free','smooth'};
        
    case 'tdalab'
        str='RIKEN TDALAB for Tensor Analysis';
        
    case 'marker'
        str='sd>op+*xv.<^ph';
        
    otherwise
        error('Unsuported input parameter.');
end
end

