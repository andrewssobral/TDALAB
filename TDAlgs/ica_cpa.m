function [S2,A2,B2] = ica_cpa(tensor,R,alfa)
% function [S2,A2,B2] = ica_cpa_v(tensor,R,alfa,truth)
%function [S2,A2,B2] = ica_cp_final2 (tensor,R,alfa)
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


    [I,J,K]=size(tensor);
    T=I*J;
    matrix=reshape(tensor, T,K);
   
    %estimate statistics
    [M,Cov]=est_cumulant(matrix,R);
for k=1:length(alfa)  
    M2=reshape(M,T,T,T);
    M2=M2(:,:,1:R);
    M2=cos(alfa(k))*M2/norm(M2(:));
    Cov=Cov/norm(Cov(:));
    M2(:,:,R+1)=sin(alfa(k))*Cov;
    %reshape to fifth order tensor
    cumulant=reshape(M2,I,J,I,J,size(M2,3));
    %estimate the decomposition with symmetry constraints
    [A{k},B{k},C{k},niter{k},Fit(k)]=PARAFAC5_ALSLS_Herm_2pairs(cumulant,R,'ALS');
    %reconstruct the mixing matrix
    Mix=kat_rao(B{k},A{k});
    %estimate the statistically independent sources
    S{k}=pinv(Mix)*matrix;
   
    clear M2
end

% Fit
% 
% for k=1:length(alfa)
%    est_ord_ampli2(S{k}.',truth) 
% end
[val,index]=min(Fit);
A2=A{index};
B2=B{index};
S2=S{index};



%**************************************************************************************************************
%----------------------------------- END OF MAIN FUNCTION ----------------------------
%**************************************************************************************************************



function [M,R] =  est_cumulant(X,nb)

% function [M] = est_cumulant(X,alfa);
%
% This function estimates the cumulant and the covariance matrix
% 


[n,T]	= size(X);

%%  source detection not implemented yet !
if nargin==1, m=n ; end;
m=n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A few parameters that could be adjusted
nem=m;

Y=X;

R	= (Y*Y' )/T ;
C	= (Y*Y.')/T ;

Yl	= zeros(1,T);
Ykl	= zeros(1,T);
Yjkl	= zeros(1,T);

Q	= zeros(m*m*m*m,1) ;
index	= 1;

for lx = 1:m ; Yl 	= Y(lx,:);
for kx = 1:m ; Ykl 	= Yl.*conj(Y(kx,:));
for jx = 1:m ; Yjkl	= Ykl.*conj(Y(jx,:));
for ix = 1:m ; 
	Q(index) = ...
	(Yjkl * Y(ix,:).')/T -  R(ix,jx)*R(lx,kx) -  R(ix,kx)*R(lx,jx) -  C(ix,lx)*conj(C(jx,kx))  ;
	index	= index + 1 ;
end ;
end ;
end ;
end




 [U,D]	= eig(reshape(Q,m*m,m*m)); 
 [la,K]	= sort(abs(diag(D)));
 
% %% reshaping the most (there are `nem' of them) significant eigenmatrice
 M	= zeros(m,nem*m);	% array to hold the significant eigen-matrices
 Z	= zeros(m)	; % buffer
 h	= m*m;
 for u=1:m:nem*m, 
 	Z(:) 		= U(:,K(h));
 	M(:,u:u+m-1)	= la(h)*Z;
 	h		= h-1; 
 end;







function [A,B,C,niter,Fit]=PARAFAC5_ALSLS_Herm_2pairs(X,R,method,tol,Niter,Ninit,A,B,C,stop_crit)
%function [A,B,C,niter,Fit]=PARAFAC5_ALSLS_Herm_2pairs(X,R,method,tol,Niter,Ninit,A,B,C,stop_crit)

%-------------------------------------------------------------------------------------------------------- 
% CANDECOMP/PARAFAC BY ALTERNATING LEAST SQUARES ALGORITHM COUPLED WITH LINE SEARCH 
% FIVE-ORDER CASE with HERMITIAN SYMMETRY on 2 pairs of modes: 
% Mode 1: A
% Mode 2: B
% Mode 3: conj(A)
% Mode 4: conj(B)
% Mode 5: C
% X is an (I1xI2xI3xI4xI5) tensor with I1=I3 and I2=I4
% Unknowns: A(I1xR), B(I2xR), C(I5xR)
%--------------------------------------------------------------------------------------------------------
% @Copyright 2009    
% Dimitri Nion 
% dimitri.nion@gmail.com
% K.U. Leuven, campus Kortrijk, Belgium
% This M-file and the code in it belongs to the holder of the copyrights. 
%-------------------------------------------------------------------------------------------------------
% Description:
% This function computes the PARAFAC decomposition of the 5th-order tensor X by means of the 
% Alternating Least Squares algorithm (ALS) possibly combined with a Line Search step to speed up convergence.
% Several Line Search techniques are possible (see references below)
% Symmetry is explicitely taken into account in ALS loops as well as in ELS stage.
%--------------------------------------------------------------------------------------------------------
% NOTE: no compression pre-processing stage is implemented in this algorithm, so the model is fitted on the 
% raw data in original space. To combine compression with symmetry, one needs specific algorithms
% that computes the best lower rank-(R1,R2,R3,R4,R5) with preservation of symmetry. See work
% of Mariya Ishteva and Lieven De Lathauwer.
% An alternative would be to ignore symmetry in the compression stage, then fit the model in compressed space
% with ALS+ELS while ignoring symmetry (so 5 different loading matrices), and finally perform a few steps in
% original space after decompression but this time with imposing symmetry.
% Also, if only the 5th dimension is long, then one could perform compression in this mode only, there is
% no problem with symmetry.
%--------------------------------------------------------------------------------------------------------
% INPUTS: 
% X: tensor of size (I1xI2xI3xI4xI5) with I1=I3 and I2=I4
% R: number of rank-1 components
%--------------------------------------------------------------------------------------------------------
% OPTIONAL INPUTS:
% - method (default 'ELS'):  string to choose the Line Search method  
%                              'ALS' (standard ALS Algorithm, i.e., no Line Search is performed)
%                              'LSH' (ALS + Line Search proposed by R.A. Harshman, i.e., STEP-size set to 1.25) 
%                              'LSB' (ALS + Line Search proposed by R. Bro, i.e. STEP-size set to n^1/3, wher n is the iteration index)
%                              'ELS' (ALS + Enhanced Line Search)
% - tol    (default 1e-6) : threshold value to stop the ALS+LS algorithm 
% - Niter  (default 5000) : max number of iterations of ALS+LS algorithm 
% - Ninit  (default 5)    : number of different random initializations tried (outputs of the algorithms are the ones stored for the best initialization)
% - A ,B ,C               : if one or more given, matrices used to initialize the algorithm and in this case, 
%                           only this initialization is tried (Ninit is set to 1)
% - stop_crit (default 'VarC') : string to choose how the algorithm is stopped
%                              'Fit' to stop when the relative decrease in Fit between 2 iterations drops below tol
%                              'VarC' to stop when the relative change of matrix C drops below tol
%------------------------------------------------------------------------------------------------------------------
% OUTPUTS:
% - A (I1xR), B(I2xR), C(I5xR): Loading matrices estimates
% - niter  : number of iterations for the ALS+LS algorithm to converge (at most Niter)
% - Fit    : Frobenius Norm of the residual tensor
% NOTE : all outputs are the ones obtained for the BEST initialization only. 
%-----------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------
% Check inputs and define default parameters
%---------------------------------------------------------------------
if length(size(X))~=5
     error('The input argument X has to be a 5th-order tensor')
end

I1=size(X,1);           % We must have size(X,1)=size(X,3)=I1 and size(X,2)=size(X,4)=I2 for the program to work
I2=size(X,2);
I3=size(X,3);
I4=size(X,4);
I5=size(X,5);
X_input=X;
if (I1~=I3)
    error('Dimensions 1 and 3 should be equal')
end
if (I2~=I4)
    error('Dimensions 2 and 4 should be equal')
end

if nargin<2
    error('At least 2 input arguments are required')
end
if exist('method')~=1 | isempty(method)
    method='ALS';
end
if exist('tol')~=1 | isempty(tol)
    tol=1e-6;
end
if exist('Niter')~=1 | isempty(Niter)
    Niter=5000;
end
if exist('Ninit')~=1 | isempty(Ninit)
    Ninit=5;
end
if exist('stop_crit')~=1 | isempty(stop_crit)
    stop_crit='VarC';
end
if isreal(X)
    data_type=1;   % we work on real data, useful to know for method ELS only
else
    data_type=2;   % work on complex data, useful to know for method ELS only
end

%-----------------------------------------------------------------------
%-------- Deal with initialization -------------------------------------
%-----------------------------------------------------------------------
init_vec=[0 0 0];  % [0 0 0] if no input matrices is given to initialize
init_given=0;
if exist('A')==1 & isempty(A)~=1
    if (size(A,1)==size(X,1)) & (size(A,2)==R)    % else a random init with the right dimensions will be built
    init_vec(1)=1; 
    A_init=A;         
    else
    disp('The input matrix A does not have the right dimensions, a random matrix will be used instead')
    end
end
if exist('B')==1 & isempty(B)~=1
    if (size(B,1)==size(X,2)) & (size(B,2)==R)
    init_vec(2)=1;
    B_init=B;
    else
    disp('The input matrix B does not have the right dimensions, a random matrix will be used instead')
    end
end
if exist('C')==1 & isempty(C)~=1
    if (size(C,1)==size(X,5)) & (size(C,2)==R)
    init_vec(3)=1;
    C_init=C;
    else
    disp('The input matrix C does not have the right dimensions, a random matrix will be used instead')
    end
end
if sum(init_vec)~=0  % at least one matrix is given with the right dimensions
    init_given=1;
    max_init=1;   
end


%---------------------------------------------------------------------
% COMPUTE THE DECOMPOSITION
%---------------------------------------------------------------------
 % Matrices that are fixed and are useful for method='ELS' only:
   Pa2=fliplr(pascal(6));
   Pa=zeros(6);
   for n=1:6
      Pa(n,n:end)=diag(Pa2,n-1).';
   end
   Ja=repmat([j^5 j^4 j^3,j^2,j,1],6,1);
   
   Pb2=fliplr(pascal(6));
   Pb=zeros(6);
   for n=1:6
      Pb(n,n:end)=diag(Pb2,n-1).';
   end 
   Jb=toeplitz([1,0,0,0,0,0],[1,j,j^2,j^3,j^4,j^5]);
   
   % Matrix unfoldings of X
    X1= reshape(permute(X,[5,4,3,2,1]),I2*I3*I4*I5,I1);  % X1 is I2I3I4I5xI1
    X2= reshape(permute(X,[1,5,4,3,2]),I3*I4*I5*I1,I2);  % X2 is I3I4I5I1xI2                   
    X3= reshape(permute(X,[2,1,5,4,3]),I4*I5*I1*I2,I3);  % X3 is I4I5I1I2xI3 
    X4= reshape(permute(X,[3,2,1,5,4]),I5*I1*I2*I3,I4);  % X4 is I5I1I2I3xI4
    X5= reshape(permute(X,[4,3,2,1,5]),I1*I2*I3*I4,I5);  % X5 is I1I2I3I4xI5
 
    Fit_best=inf;              % Useful to select the best initialization

    for ninit=1:Ninit   % Try Ninit different initializations
            
            %--------------------------------------------------------------------------------------------------
            % Build INITIALIZATION
            %---------------------------------------------------------------------------------------------------
            if init_given==0     
                [A,B,C]=PARAFAC_init(I1,I2,I5,R,data_type);
            elseif init_given==1
                    if sum(init_vec)==1     % only one matrix was provided
                         [A,B,C]=PARAFAC_init(I1,I2,I5,R,data_type);  % generate random initializations and then overwrite by the one given in input
                         if init_vec(1)==1   % only A was provided
                             A = A_init;
                             B = X2.'/kat_rao(A,kat_rao(B,kat_rao(C,A))).';  % use X2  I3I4I5I1xI2
                             C = X5.'/kat_rao(A,kat_rao(B,kat_rao(A,B))).';  % use X5 I1I2I3I4xI5
                         elseif init_vec(2)==1  % only B was provided
                            B = B_init;
                            A = X1.'/kat_rao(B,kat_rao(A,kat_rao(B,C))).';  % use X1  I2I3I4I5xI1
                            C = X5.'/kat_rao(A,kat_rao(B,kat_rao(A,B))).';  % use X5 I1I2I3I4xI5
                         elseif init_vec(3)==1   % only C was provided
                            C = C_init;
                            A = X1.'/kat_rao(B,kat_rao(A,kat_rao(B,C))).';
                            B = X2.'/kat_rao(A,kat_rao(B,kat_rao(C,A))).';
                         end

                    elseif sum(init_vec)==2  % two matrices were provided
                        [A,B,C]=PARAFAC_init(I1,I2,I5,R,data_type);  % generate random initializations and then overwrite by the one given in input
                         if init_vec(1)==0    % B and C are given and A=[]
                           B = B_init;C = C_init;
                           A = X1.'/kat_rao(B,kat_rao(A,kat_rao(B,C))).';
                         elseif init_vec(2)==0      % A and C are given in input and B=[]
                           A = A_init; C = C_init;
                           B = X2.'/kat_rao(A,kat_rao(B,kat_rao(C,A))).';   % so we update B to exploit A and C
                         elseif init_vec(3)==0  % A and B were given 
                           A = A_init; B = B_init; 
                           C = X5.'/kat_rao(A,kat_rao(B,kat_rao(A,B))).';  % so we update C  to exploit A and B 
                         end
                    elseif sum(init_vec)==3  % three matrices were provided
                          A = A_init; B = B_init; C = C_init;
                    end
            end
            
                % init useful for Line Search
                A1=A;A2=A;
                B1=B;B2=B;
                C1=C;C2=C;  

            
            %-------------------------------------------------------------------------
            % LOOP for alternating updates
            %-------------------------------------------------------------------------
            Fit=norm(X5-kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B))))*C.','fro');
           
            VarC = norm(C*diag(1./sqrt(sum(C.*conj(C)))));
            stop=0;
            Fit_vec=[];
            niter=0;
   
            while stop==0
        
                niter=niter+1;
                VarC_old=VarC;
                Fit_old=Fit;
                
               %----------------------
               % Perform Line Search
               %----------------------
                  [A,B,C] = Line_Search(A1,B1,C1,A2,B2,C2,X,method,niter,data_type,Pa,Ja,Pb,Jb); 
             
               %--------------------------------
               % Perform ALS updates
               %--------------------------------- 
             
                    A = (inv((B'*B).*(conj(B')*conj(B)).*(conj(A')*conj(A)).*(C'*C)) * (kat_rao(kat_rao(B,conj(A)),kat_rao(conj(B),C))'*X1) ).';
                    B = (inv((A'*A).*(conj(A')*conj(A)).*(conj(B)'*conj(B)).*(C'*C)) * (kat_rao(kat_rao(conj(A),conj(B)),kat_rao(C,A))'*X2) ).';
                    C = (inv((A'*A).*(conj(A')*conj(A)).*(B'*B).*(conj(B)'*conj(B))) * (kat_rao(kat_rao(A,B),kat_rao(conj(A),conj(B)))'*X5) ).';
                    A = (inv((B'*B).*(conj(B')*conj(B)).*(A'*A).*(C'*C)) * (kat_rao(kat_rao(conj(B),C),kat_rao(A,B))'*X3) ).';
                    A = conj(A);
                    B = (inv((A'*A).*(conj(A')*conj(A)).*(B'*B).*(C'*C)) * (kat_rao(kat_rao(C,A),kat_rao(B,conj(A)))'*X4) ).';
                    B = conj(B);
                    
%----- By exploiting properties of khatri-rao product, the latter implementation is much faster than             
%                          A = X1.'/kat_rao(B,kat_rao(conj(A),kat_rao(conj(B),C))).';  % use X1  I2I3I4I5xI1
%                          B = X2.'/kat_rao(conj(A),kat_rao(conj(B),kat_rao(C,A))).';  % use X2  I3I4I5I1xI2
%                          C = X5.'/kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B)))).';  % use X5 I1I2I3I4xI5
%                          A = X3.'/kat_rao(conj(B),kat_rao(C,kat_rao(A,B))).';  % use X3 I4I5I1I2xI3
%                          A=conj(A);
%                          B = X4.'/kat_rao(C,kat_rao(A,kat_rao(B,conj(A)))).';  % use X4 I5I1I2I3xI4
%                          B=conj(B);

               %----------------------------------
               % Normalization to avoid overflow
               %----------------------------------
                        normA=sqrt(sum(A.*conj(A)));  % gives the Frobenius norm of each column of A
                        normB=sqrt(sum(B.*conj(B)));
                        normC=sqrt(sum(C.*conj(C)));
                        prod_norm=normA.*normA.*normB.*normB.*normC;
                        Scale_mat=diag(prod_norm.^(1/5)); 
                        % to get an equal repartition of the power of each rank-1 tensor over the 3 vectors:
                        A=A*diag(1./normA)*Scale_mat;
                        B=B*diag(1./normB)*Scale_mat;
                        C=C*diag(1./normC)*Scale_mat;
                        
                %--------------------------------------------------------------
                % STOP CRITERION
                %---------------------------------------------------------------
                if strcmp(stop_crit,'Fit')==1
                    Fit=norm(X5-kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B))))*C.','fro');
                       if (abs((Fit-Fit_old)/Fit_old)<tol) | (niter==Niter) | (Fit<tol)
                           stop=1;
                       end
                elseif strcmp(stop_crit,'VarC')==1  
                    VarC = norm(C*diag(1./sqrt(sum(C.*conj(C))))-C1*diag(1./sqrt(sum(C1.*conj(C1)))));
                       if (abs((VarC-VarC_old)/VarC_old)<tol) | (niter==Niter) | (VarC<tol)
                           stop=1;
                       end
                end
                       
                %-------------------------------------------------       
                % Store matrices to prepare next Line Search step
                %-------------------------------------------------
                        A2=A1;B2=B1;C2=C1;
                        A1=A;B1=B;C1=C;
                        
                       
                                               
            end    % Algorithm has converged for this initialization
            
           % Compute the Fit
            Fit=norm(X5-kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B))))*C.','fro');
    
          % Select this init if it is better than previous one
                     if (Fit < Fit_best)
                        Fit_best=Fit;
                        A_best=A;
                        B_best=B;
                        C_best=C;
                        niter_best=niter;
                        Fit_best=Fit;
                     end
     end   % END of ALL init
     
     % Output arguments corresponding to the best initialization
            A=A_best;
            B=B_best;
            C=C_best;
            niter=niter_best;
            Fit=Fit_best;
            
            
            

        
%*********************************************************************************************************        
    function [A,B,C]=PARAFAC_init(I1,I2,I5,R,data_type);
               if data_type==1      % real initialization
                    A=randn(I1,R);
                    B=randn(I2,R);
                    C=randn(I5,R);
                    if I1>=R;A=orth(A);else;A=orth(A')';end
                    if I2>=R;B=orth(B);else;B=orth(B')';end
                    if I5>=R;C=orth(C);else;C=orth(C')';end
                elseif data_type==2  % complex init
                    A=randn(I1,R)+j*randn(I1,R);
                    B=randn(I2,R)+j*randn(I2,R);
                    C=randn(I5,R)+j*randn(I5,R);
                    if I1>=R;A=orth(A);else;A=orth(A')';end
                    if I2>=R;B=orth(B);else;B=orth(B')';end
                    if I5>=R;C=orth(C);else;C=orth(C')';end
                end
        
%**************************************************************************************************        
function C = kat_rao(A,B)
%--------------------------------------------------------------------------------------------------
% Calculate the Kathri Rao product of 2 matrices
%--------------------------------------------------------------------------------------------------
        I=size(A,1);R1=size(A,2);
        J=size(B,1);R2=size(B,2);
        if R1~=R2
            error('Input matrices must have the same number of columns for Khatri Rao product')
        end
        
        C=zeros(I*J,R1);
        for j=1:R1
            C(:,j)=reshape(B(:,j)*A(:,j).',I*J,1);
        end
 

%************************************************************************************************************        
% START of functions useful only for Line Search
%************************************************************************************************************
function [A,B,C] = Line_Search(A1,B1,C1,A2,B2,C2,X,method,niter,data_type,Pa,Ja,Pb,Jb) 
        %-----------------------------------------------------------------------------------------------------
        % This function performs Line Search for PARAFAC, Five order case with symmetry on 2 pairs of modes
        % Mode 1: A
        % Mode 2: B
        % Mode 3: conj(A)
        % Mode 4: conj(B)
        % Mode 5: C
        %-----------------------------------------------------------------------------------------------------
        % @Copyright Dimitri Nion
        % Version: May 2009
        %-----------------------------------------------------------------------------------------------------
        % INPUTS: - Old estimates of the loading matrices (A1,A2,B1,B2,C1,C2)
        %         - X id the observed tensor X
        %         - method = 'ALS', 'LSH','LSB' or 'ELS'
        %         - niter is the current ALS iteration (useful for LSB only)
        %          
        %         - data_type = 1 for real data or 2 for complex, useful for ELS method only
        %------------------------------------------------------------------------------------------------------
        % OUTPUTS: Interpolated loading matrices ready to be used in the next ALS update
        %------------------------------------------------------------------------------------------------------
        % Line Search at the n th ALS iteration consists of interpolation of the loading 
        % matrices A, B, C from their estimates at (n-1)th iteration (given by A1,B1,C1) 
        % and at (n-2)th iteration (given by A2,B2,C2).
        % The crucial point is to find a good step size in the interpolation directions
        % This step size is denoted STEP here
        %-------------------------------------------------------------------------------------------------------
        % If method = 'ALS' then no interpolation is done, i.e. we have the
        % standard ALS algorithm. This corresponds to a step-size STEP=1
        %-------------------------------------------------------------------------------------------------------
        % If method = 'LSH' 
        % then we perform Line Search proposed by R.A. Harshman. STEP is fixed to 1.25 and 
        % if the interpolated matrices did not allow a decrease in Fit, then we annihilate 
        % interpolation by choosing STEP=1 (i.e. a standard ALS update is performed)
        % [5]: R.A. Harshman, "Foundations of the PARAFAC procedure: Model and Conditions for
        %      an explanatory Multi-mode Factor Analysis", UCLA Working Papers in
        %      Phonetics, vol.16, pp 1-84, 1970
        %-------------------------------------------------------------------------------------------------------
        % If method = 'LSB' 
        % then we perform Line Search proposed R. Bro. STEP is fixed to niter^(1/3), where niter 
        % is the current ALS iteration and then we validate interpolation only if we got a decrease 
        % in Fit, otherwise we annihilate Line Search by choosing STEP=1
        % [2] R. Bro, "Multi-Way Analysis in the Food Industry: Models, Algorithms, and
        %     Applications", PhD. dissertation, University of Amsterdam, 1998
        %--------------------------------------------------------------------------------------------------------
        % If method = 'ELS'
        % then we perform Enhanced Line Search which looks for the STEP value that gives the best 
        % decrease in Fit, i.e., the cost function is regarded as a function depending on STEP only,
        % where the loading matrices from the two previous iterations are fixed.
        % There are 2 cases:
        % - For REAL valued tensors, we look for the best REAL STEP, as explained in [6-7], where Enhanced 
        %   Line Search was introduced.
        % - For COMPLEX tensors, we look for the best COMPLEX STEP, as proposed in [8], where the work of [6-7]
        %   was extended to the complex case.
        % Note that, in principle, the most efficient interpolation would be to compute a different
        % STEP for each of the three interpolation directions (see [6-7] for details). Here, we consider the same STEP
        % in each direction, which gives satisfaction in practice. Indeed, the purpose is just to speed up the convergence
        % of the standard ALS with a Line Search method that has a low complexity.
        %
        % [7]  M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A Novel Method to Accelerate 
        %      PARAFAC", SIAM Journal Matrix Anal. and Appl. (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008 
        % [8]: D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for Complex-Valued Tensor Decompositions. 
        %      Application in DS-CDMA", Signal Processing, vol. 88, issue 3, pp. 749-755, March 2008.
        %----------------------------------------------------------------------------------------------------------------
        niter_start=4;   % Line Search will start after niter_start iterations (at least 2)

        % First compute the interpolation directions
        Ga=A1-A2;
        Gb=B1-B2;
        Gc=C1-C2;
        
        I1=size(X,1);           % We must have size(X,1)=size(X,3)=I1 and size(X,2)=size(X,4)=I2 for the program to work
        I2=size(X,2);
        I3=size(X,3);
        I4=size(X,4);
        I5=size(X,5);
        R=size(A1,2);
        X_mat= reshape(permute(X,[4,3,2,1,5]),I1*I2*I3*I4,I5);  % X5 is I1I2I3I4xI5
        
        %-----------------------------------------------------
        % ALS             
        %-----------------------------------------------------
        if   strcmp(method,'ALS')==1
             A=A1;B=B1;C=C1;     % Matrices are not interpolated

        %-----------------------------------------------------
        % LSH             
        %-----------------------------------------------------
        elseif strcmp(method,'LSH')==1
                if niter<niter_start
                    A=A1;B=B1;C=C1;
                else
                    % Compute the current Fit to the model
                    Fit=norm(X_mat-kat_rao(A1,kat_rao(B1,kat_rao(conj(A1),conj(B1))))*C1.','fro');
                   
                    % Interpolate
                        STEP=1.25;
                        A=A2+STEP*Ga;
                        B=B2+STEP*Gb;
                        C=C2+STEP*Gc;
                    % Compute Fit after Interpolation
                    Fit_new=norm(X_mat-kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B))))*C.','fro');
                    % Accept or Reject Interpolation                   
                        if Fit_new > Fit     % reject  
                            A=A1;B=B1;C=C1;
                        end
                end   
        %-----------------------------------------------------
        % LSB          
        %-----------------------------------------------------
        elseif strcmp(method,'LSB')==1
                if niter<niter_start
                    A=A1;B=B1;C=C1;
                else
                    % Compute the current Fit to the model
                        Fit=norm(X_mat-kat_rao(A1,kat_rao(B1,kat_rao(conj(A1),conj(B1))))*C1.','fro');
                    % Interpolate
                        STEP=niter^(1/3);
                        A=A2+STEP*Ga;
                        B=B2+STEP*Gb;
                        C=C2+STEP*Gc;
                    % Compute Fit after Interpolation
                    Fit_new=norm(X_mat-kat_rao(A,kat_rao(B,kat_rao(conj(A),conj(B))))*C.','fro');
                    % Accept or Reject Interpolation                   
                        if Fit_new > Fit     % reject  
                            A=A1;B=B1;C=C1;
                        end
                end

        %-----------------------------------------------------
        % ELS     2 cases: real or complex data
        %-----------------------------------------------------
        elseif strcmp(method,'ELS')==1

                if niter<niter_start
                    A=A1;B=B1;C=C1;     
                else
                   % Interpolate
                        STEP = Find_STEP_ELS(A2,B2,C2,Ga,Gb,Gc,X_mat,data_type,Pa,Ja,Pb,Jb);   
                        A=A2+STEP*Ga;
                        B=B2+STEP*Gb;
                        C=C2+STEP*Gc;
                                         
                end
        end

        % END OF Line Search Function
        %**********************************************************************

%****************************************************************************************************************
function [STEP]=Find_STEP_ELS(A2,B2,C2,Ga,Gb,Gc,X_mat,data_type,Pa,Ja,Pb,Jb,tol,max_iter_in)
% This is the core function of Enhanced Line Search. It computes STEP in the case
% of 5-way PARAFAC with symmetry on 2 pairs of modes
% Mode 1: A
% Mode 2: B
% Mode 3: conj(A)
% Mode 4: conj(B)
% Mode 5: C
%-----------------------------------------------------------------------------
% [7]  M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A Novel Method to Accelerate 
%     PARAFAC", SIAM Journal Matrix Anal. and Appl. (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008 
% [8]: D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for Complex-Valued Tensor Decompositions. 
%      Application in DS-CDMA", Signal Processing, vol. 88, issue 3, pp. 749-755, March 2008.
%-------------------------------------------------------------------------------------------------------

if exist('tol')~=1 | isempty(tol)    % tolerance for the inner loop of ELS with complex step
    tol=1e-4;
end 

if exist('max_iter_in')~=1 | isempty(max_iter_in)    % max nb of iterations for the inner loop of ELS with complex step
    max_iter_in=10;
end 
   
    
          Ga_Gb=kat_rao(Ga,Gb);
          a_b=kat_rao(A2,B2);
          Ga_b=kat_rao(Ga,B2);
          a_Gb=kat_rao(A2,Gb);
          %---- 4 G terms
          Ga_Gb_Ga_Gb = kat_rao(Ga_Gb,conj(Ga_Gb));
          %----- 3 G terms
          Ga_Gb_Ga_b  = kat_rao(Ga_Gb,conj(Ga_b));
          Ga_Gb_a_Gb  = kat_rao(Ga_Gb,conj(a_Gb));
          Ga_b_Ga_Gb  = kat_rao(Ga_b,conj(Ga_Gb));
          a_Gb_Ga_Gb  = kat_rao(a_Gb,conj(Ga_Gb));
          %------ 2 G terms
          Ga_Gb_a_b  = kat_rao(Ga_Gb,conj(a_b));
          Ga_b_Ga_b  = kat_rao(Ga_b,conj(Ga_b));
          Ga_b_a_Gb  = kat_rao(Ga_b,conj(a_Gb));
          a_Gb_Ga_b  = kat_rao(a_Gb,conj(Ga_b));
          a_Gb_a_Gb  = kat_rao(a_Gb,conj(a_Gb));
          a_b_Ga_Gb  = kat_rao(a_b,conj(Ga_Gb));
          %-------1 G term
          Ga_b_a_b  = kat_rao(Ga_b,conj(a_b));
          a_Gb_a_b  = kat_rao(a_Gb,conj(a_b));
          a_b_Ga_b  = kat_rao(a_b,conj(Ga_b));
          a_b_a_Gb  = kat_rao(a_b,conj(a_Gb));
          %----- 0 G term
          a_b_a_b   = kat_rao(a_b,conj(a_b));
          
          Mat5 = Ga_Gb_Ga_Gb*Gc.';      %---> coeff of STEP^5 
          Mat4 = Ga_Gb_Ga_Gb*C2.' ...   %---> coeff of STEP^4
               + (a_Gb_Ga_Gb ...
               + Ga_b_Ga_Gb ...
               + Ga_Gb_a_Gb ...
               + Ga_Gb_Ga_b)*Gc.';
          Mat3 = (Ga_Gb_Ga_b ...     %---> coeff of STEP^3
               + Ga_Gb_a_Gb ...
               + Ga_b_Ga_Gb ...
               + a_Gb_Ga_Gb)*C2.' ...
               + (Ga_Gb_a_b ...
               + Ga_b_Ga_b ...
               + Ga_b_a_Gb ...
               + a_Gb_Ga_b ...
               + a_Gb_a_Gb ...
               + a_b_Ga_Gb)*Gc.';
          Mat2 = (a_b_a_Gb ...     %---> coeff of STEP^2
               + a_b_Ga_b ...
               + Ga_b_a_b ...
               + a_Gb_a_b)*Gc.' ...
               + (a_Gb_a_Gb ...
               + a_Gb_Ga_b ...
               + a_b_Ga_Gb ...    
               + Ga_b_a_Gb ...
               + Ga_b_Ga_b ...
               + Ga_Gb_a_b)*C2.';
          Mat1 = a_b_a_b*Gc.' ...   %---> coeff of STEP
               + (Ga_b_a_b ...
               + a_Gb_a_b ...
               + a_b_Ga_b ... 
               + a_b_a_Gb)*C2.';
          Mat0 = a_b_a_b*C2.'-X_mat;      %---> constant term
          
          M=[Mat5(:) Mat4(:) Mat3(:) Mat2(:) Mat1(:) Mat0(:)]; 

        %--------------------------REAL CASE-----------------------------------------
        % Method by M. Rajih and P. Comon
                       if data_type==1  % choose the best real STEP (Rajih and Comon method)

                               H_mat=M'*M;  
                               % Now we define the coefficients of the 10th order polynomial, which correspond to the sum of the different
                               % antidiagonal terms of H_mat
                               H_mat2=fliplr(H_mat);  % anti-diagonal of H_mat is now diagonal of H_mat2
                               k=0;
                               for n=5:-1:-5
                               k=k+1;
                               pol(1,k)=sum(diag(H_mat2,n));
                               end
                               pol_der=polyder(pol);
                               sqrts=roots(pol_der);

                               %Take modulus (the step size has to be real)                           
                                sqrts=[abs(sqrts);1];   % 1 is added to the set of possible steps (it corresponds to ALS without line search)                          
                               % Choice of R_opt
                                   extremum=polyval(pol,abs(sqrts));
                                   STEP=sqrts(find(extremum==min(extremum),1));
                               
                               

        %--------------------------COMPLEX CASE----------------------------------------                             
        % Method by D. Nion and L. De Lathauwer
                       elseif data_type==2  
                            H_mat=M'*M;  
                           
                           %--------------
                            % Init
                            %-------------
                            STEP=1;
                            a=1;  % real part of step
                            b=0;  % imag part of step
                            
                            u=[STEP^5 STEP^4 STEP^3 STEP^2 STEP 1].';
                            norm_new=abs(u'*H_mat*u);
                            norm_diff=norm_new;
                            niter_in=0;

                            % Alternate update of real part a and imag part b of STEP 
                             while (norm_diff > tol) & (niter_in < max_iter_in)
                                    niter_in=niter_in+1;
                                    norm_old=norm_new;
                                %---------------------
                                % Update imag part b  
                                %---------------------
                                 [Da] = Build_Da(a,Pa,Ja);
                                 Ha_mat=Da'*H_mat*Da;
                                 
                                 % Now we define the coefficients of the 6th order polynomial
                                 Ha_mat2=fliplr(Ha_mat);  % anti-diagonal of H_mat is now diagonal of H_mat2
                                   k=0;
                                   for n=5:-1:-5
                                   k=k+1;
                                   pol(1,k)=sum(diag(Ha_mat2,n));
                                   end
                                pol_der=polyder(pol);
                                sqrts=roots(pol_der);
                               
                               sqrts2=sqrts;
                               for nn=1:length(sqrts) 
                                if imag(sqrts(nn))~=0
                                   sqrts(nn)=real(sqrts(nn));   % replace by real part
                                   sqrts2(nn)=abs(sqrts(nn));   % replace by modulus
                                end
                               end
                               sqrts=[sqrts;sqrts2];  % we add 1 (i.e. annihilate Line Search) in the set of possible STEPS values
                               % Choice of optimal STEP
                                   extremum=polyval(pol,sqrts);
                                   b=sqrts(find(extremum==min(extremum),1));
                                   b=b(1);
                               
                               %---------------------
                                % Update real part a  
                                %---------------------
                                 [Db] = Build_Db(b,Pb,Jb);
                                 Hb_mat=Db'*H_mat*Db;
                                 
                                 % Now we define the coefficients of the 6th order polynomial
                                 Hb_mat2=fliplr(Hb_mat);  % anti-diagonal of H_mat is now diagonal of H_mat2
                                   k=0;  
                                   for n=5:-1:-5
                                   k=k+1;
                                   pol(1,k)=sum(diag(Hb_mat2,n));
                                   end
                                pol_der=polyder(pol);
                                sqrts=roots(pol_der);
                               
                               sqrts2=sqrts;
                               for nn=1:length(sqrts) 
                                if imag(sqrts(nn))~=0
                                   sqrts(nn)=real(sqrts(nn));   % replace by real part
                                   sqrts2(nn)=abs(sqrts(nn));   % replace by modulus
                                end
                               end
                               sqrts=[sqrts;sqrts2];  % we add 1 (i.e. annihilate Line Search) in the set of possible STEPS values
                               % Choice of optimal STEP
                                   extremum=polyval(pol,sqrts);
                                   a=sqrts(find(extremum==min(extremum),1));
                                   a=a(1);


                                 %***********************************************
                                 %updated STEP
                                 STEP=a+j*b;
                                 u=[STEP^5 STEP^4 STEP^3 STEP^2 STEP 1].';
                                 norm_new=abs(u'*H_mat*u);
                                 norm_diff=abs(norm_new-norm_old);
                             end 
                       end
                       
%-------------------------------------------------------------------------------------------------
    function [Da] = Build_Da(a,Pa,Ja)
        % Given that u =[STEP^5 STEP^4 STEP^3 STEP^2 STEP 1].' and that STEP is complex valued: STEP=a+ib
        % This function builds the 6x6 matrix Da such that u = Da*[b^5 b^4 b^3 b^2 b 1].' , where Da is built from a only.
        % Indeed Da is upper triangular and its coefficients are obtained from pascal triangle
        a_mat=toeplitz([1,0,0,0,0,0],[1 a a^2 a^3 a^4 a^5]);
        Da=Pa.*Ja.*a_mat;
        
%--------------------------------------------------------------------------------------------        
         function [Db] = Build_Db(b,Pb,Jb)
        % Given that u =[STEP^5 STEP^4 STEP^3 STEP^2 STEP 1].' and that STEP is complex valued: STEP=a+ib
        % This function builds the 6x6 matrix Db such that u = Db*[a^5 a^4 a^3 a^2 a 1].' , where Db is built from b only.
        % Indeed Db is upper triangular and its coefficients are given by pascal triangle
        b_mat=toeplitz([1,0,0,0,0,0],[1 b b^2 b^3 b^4 b^5]);
        Db=Pb.*Jb.*b_mat;        

%**************************************************************************        
% END of functions useful only for ALS and ALS+Line Search
%**************************************************************************                     
m b only.
        % Indeed Db is upper triangular and its coefficients are given by pascal triangle
        b_mat=toeplitz([1,0,0,0,0,0],[1 b b^2 b^3 b^4 b^5]);
        Db=Pb.*Jb.*b_mat;        

%**************************************************************************        
% END of functions useful only for ALS and ALS+Line Search
%**************************************************************************                     



