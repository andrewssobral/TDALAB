function [Y,A,phi_vec]=cp4_alsls(X,R,lsearch,Tol,MaxIt,Ninit)
% function [A1,A2,A3,A4,it,phi,phi_vec]=cp4_alsls(X,R,lsearch,Tol,MaxIt,Ninit)
%CP4_ALSLS CANDECOMP/PARAFAC decomposition of a fourth-order tensor (CP4).
%   [A1,A2,A3,A4] = cp4_alsls(X,R) computes a CANDECOMP/PARAFAC decomposition
%   of a fourth-order tensor X in R rank-one terms, stored in the factor 
%   matrices A1, A2, A3, A4, belonging to the first, second, third and fourth 
%   mode, respectively.
%   The alternating least squares (ALS) algorithm is used, possibly coupled
%   with different choices of line search to speed up convergence.
%
% INPUTS: 
% X: tensor of size (I1xI2xI3xI4) 
% R: number of rank-1 components
%
% OPTIONAL INPUTS:
% - lsearch  (default 'lsb') : The type of line search to use. 
%           The different types are:
%            'none'  (standard ALS Algorithm, i.e., no Line Search is performed)
%             'lsh'  (Line Search proposed by Harshman in [5], i.e., 
%                     STEP-size set to 1.25) 
%             'lsb'  (Line Search proposed by Bro in [2], i.e. STEP-size set to 
%                     n^1/3, wher n is the iteration index)
%             'elsr' (Exact Line Search with optimal real-valued step [6,7])
%             'elsc' (Exact Line Search with optimal complex-valued step [6,7])
% - Tol    (default 1e-4) : threshold value to stop the algorithm 
% - MaxIt  (default 1000) : max number of iterations of the algorithm 
% - Ninit  (default 3)    : number of different random initializations tried 
%
% OUTPUTS:
% - A1 (I1xR), A2(I2xR), A3(I3xR), A4(I4xR): Loading matrices estimates
% - it      : number of iterations for the algorithm to converge (at most MaxIt)
% - phi     : Frobenius Norm of the residual tensor
% - phi_vec : Evolution of phi
%  NOTE: all outputs are the ones obtained for the BEST initialization only, the
%  latter being selected as the one that yields the smallest final value of phi.
%-------------------------------------------------------------------------------
% Author: Dimitri Nion   (Feedback: dimitri.nion@gmail.com)
% @Copyright May 2010
% All rights reserved. This M-file and the code in it belongs to the holder 
% of the copyright. For non commercial use only.
%-------------------------------------------------------------------------------

kr=@khatrirao;
% Check inputs and define default parameters
if length(size(X))~=4;
    error('The input argument X has to be a 4th-order tensor')
end

I1=size(X,1);          
I2=size(X,2);
I3=size(X,3);
I4=size(X,4);

if nargin<2;error('At least 2 input arguments are required');end
if exist('lsearch')~=1 || isempty(lsearch);         lsearch='lsb';          end
if exist('Tol')~=1 || isempty(Tol);                 Tol=1e-6;                end
if exist('MaxIt')~=1 || isempty(MaxIt);             MaxIt=1e-4;              end
if exist('Ninit')~=1 || isempty(Ninit);             Ninit=3;                 end


%---------------------------------------------------------------------
% COMPUTE THE DECOMPOSITION
%---------------------------------------------------------------------
phi_best=inf;              % Useful to select the best initialization

% Matrix unfoldings of X
X1= reshape(permute(X,[4,3,2,1]),I2*I3*I4,I1);  % X1 is I2I3I4xI1
X2= reshape(permute(X,[1,4,3,2]),I3*I4*I1,I2);  % X2 is I3I4I1xI2                   
X3= reshape(permute(X,[2,1,4,3]),I4*I1*I2,I3);  % X3 is I4I1I2xI3 
X4= reshape(permute(X,[3,2,1,4]),I1*I2*I3,I4);  % X4 is I1I2I3xI4

for ninit=1:Ninit   % Try max_init different initializations

    display('-------------------------------------------------------------')
    display(['-------- Running cp4_alsls - initialization ',num2str(ninit)])
    
%---------------------------------------------------------------------------
% Build INITIALIZATION
%--------------------------------------------------------------------------- 
if isreal(X); data_type='real'; else; data_type='complex';end
[A1,A2,A3,A4]=CP4_init(I1,I2,I3,I4,R,data_type);

% init useful for Line Search
A1p=A1;A1q=A1;
A2p=A2;A2q=A2;
A3p=A3;A3q=A3;
A4p=A4;A4q=A4;

%-------------------------------------------------------------------------
% LOOP for alternating updates
%-------------------------------------------------------------------------
phi=norm(X4-kr(A1,A2,A3)*A4.','fro');
stop=0;
it=0;
phi_vec=[];

while stop==0
    it=it+1;
    phi_old=phi;

   %------------------------------------------------------------------------
   % LINE SEARCH STEP: INTERPOLATE MATRICES
   %------------------------------------------------------------------------
   [A1,A2,A3,A4] = Line_Search4way(A1p,A2p,A3p,A4p,A1q,A2q,A3q,A4q,X4,phi,lsearch,it);

   %--------------------------------
   % Perform ALS updates
   %--------------------------------- 
   A1 = (inv((A2'*A2).*(A3'*A3).*(A4'*A4))* (kr(A2,A3,A4)'*X1) ).';   % X1 is I2I3I4xI1
   A2 = (inv((A3'*A3).*(A4'*A4).*(A1'*A1))* (kr(A3,A4,A1)'*X2) ).';   % X2 is I3I4I1xI2  
   A3 = (inv((A4'*A4).*(A1'*A1).*(A2'*A2))* (kr(A4,A1,A2)'*X3) ).';   % X3 is I4I1I2xI3 
   A4 = (inv((A1'*A1).*(A2'*A2).*(A3'*A3))* (kr(A1,A2,A3)'*X4) ).';   % X4 is I1I2I3xI4

   %----------------------------------
   % Normalization to avoid overflow
   %----------------------------------
    normA1=sqrt(sum(A1.*conj(A1)));  % gives the Frobenius norm of each column of A
    normA2=sqrt(sum(A2.*conj(A2)));
    normA3=sqrt(sum(A3.*conj(A3)));
    normA4=sqrt(sum(A4.*conj(A4)));
    prod_norm=normA1.*normA2.*normA3.*normA4;
    Scale_mat=diag(prod_norm.^(1/4)); 
    % to get an equal repartition of the power of each rank-1 tensor over the 5 modes:
    A1=A1*diag(1./normA1)*Scale_mat;
    A2=A2*diag(1./normA2)*Scale_mat;
    A3=A3*diag(1./normA3)*Scale_mat;
    A4=A4*diag(1./normA4)*Scale_mat;

    %--------------------------------------------------------------
    % Calculate the new phi  and decide to stop or not
    %---------------------------------------------------------------
    % ---- STOP criterion: evolution of phi -----
    phi=norm(X4-kr(A1,A2,A3)*A4.','fro');
    phi_vec=[phi_vec,phi];
    if (abs((phi-phi_old)/phi_old)<Tol) || (it==MaxIt) || (phi<Tol)
       stop=1;
    end

     %---- Other criterion: 
     % evolution of matrix A1 (to avoid computation of phi in the loop) -----
     % VarA1 = norm(A1*diag(1./sqrt(sum(A1.*conj(A1))))-A1old*diag(1./sqrt(sum(A1old.*conj(A1old)))));
     % if (abs((VarA1-VarA1_old)/VarA1_old)<Tol) || (it==MaxIt) || (VarA1<Tol)
     %         stop=1;
     % end

    %-------------------------------------------------       
    % Store matrices to prepare next Line Search step
    %-------------------------------------------------
    A1q=A1p;A1p=A1;
    A2q=A2p;A2p=A2;
    A3q=A3p;A3p=A3;
    A4q=A4p;A4p=A4;


end    % Algorithm has converged for this initialization

display(['-------- Stopped. Fit= ',num2str(phi),' Iterations= ',num2str(it)])
display('-------------------------------------------------------------') 
% Select this init if it is better than previous one
         if (phi < phi_best)
            phi_best=phi;
            phi_vec_best=phi_vec;
            A1_best=A1;
            A2_best=A2;
            A3_best=A3;
            A4_best=A4;
            it_best=it;
            phi_vec=[];
         end
end   % END of ALL init
     

% Output arguments corresponding to the best initialization
A1=A1_best;
A2=A2_best;
A3=A3_best;
A4=A4_best;
it=it_best;
phi=phi_best;
phi_vec=phi_vec_best;

Y=ktensor(ones(R,1),A1,A2,A3,A4);
A=Y.U(:);
end         
        
%*******************************************************************************        
function [A1,A2,A3,A4]=CP4_init(I1,I2,I3,I4,R,data_type)
% Initialize CP4 with orthogonal random matrices
if strcmp(data_type,'real')==1      % real initialization
    A1=randn(I1,R);
    A2=randn(I2,R);
    A3=randn(I3,R);
    A4=randn(I4,R);
    if I1>=R;A1=orth(A1);else;A1=orth(A1')';end
    if I2>=R;A2=orth(A2);else;A2=orth(A2')';end
    if I3>=R;A3=orth(A3);else;A3=orth(A3')';end
    if I4>=R;A4=orth(A4);else;A4=orth(A4')';end

elseif strcmp(data_type,'complex')==1  % complex init
    A1=randn(I1,R)+j*randn(I1,R);
    A2=randn(I2,R)+j*randn(I2,R);
    A3=randn(I3,R)+j*randn(I3,R);
    A4=randn(I4,R)+j*randn(I4,R);
    if I1>=R;A1=orth(A1);else;A1=orth(A1')';end
    if I2>=R;A2=orth(A2);else;A2=orth(A2')';end
    if I3>=R;A3=orth(A3);else;A3=orth(A3')';end
    if I4>=R;A4=orth(A4);else;A4=orth(A4')';end
end
end
               
%*******************************************************************************        
function C = kr(A1,A2,varargin)
% Kathri Rao product of N matrices
% INPUTS:
% - Matrices A1 (I1xR) and A2 (I2xR)
% OPTIONAL INPUTS:
% - varargin is a list of matrices A3,A4,A5,..., all of them having R columns
% OUTPUT:
% - Matrix C is the Khatri-product of A1 and A2 (or A1,A2,A3,A4,A5) 
N=size(varargin,2);
I1=size(A1,1);
I2=size(A2,1);
R1=size(A1,2);
R2=size(A2,2);
if R1~=R2
    error('Input matrices must have the same number of columns')
end
% perform Khathri-Rao product of A1 and A2
C=zeros(I1*I2,R1);
for r=1:R1
    C(:,r)=reshape(A2(:,r)*A1(:,r).',I1*I2,1);
end
% Do a recursive call to compute the product with the remaining matrices
for p=1:N
    C=kr(C,varargin{p});
end
end   
    
%**************************************************************************************************             
function [A1,A2,A3,A4] = Line_Search4way(A1p,A2p,A3p,A4p,A1q,A2q,A3q,A4q,X4,phi,lsearch,it)
% This function performs Line Search for 4-way PARAFAC
% Mode 1: A1
% Mode 2: A2
% Mode 3: A3
% Mode 4: A4
%
% INPUTS: - Estimates of the loading matrices, one iteration ago:
%           A1p,A2p,A3p,A4p
%         - Estimates of the loading matrices, two iterations ago:
%           A1q,A2q,A3q,A4q
%         - X4 is the I1I2I3xI4 matrix unfolding of the observed tensor X
%         - phi  is the current value of phi (i.e., Frobenius norm of residual tensor at the current iteration)
%         - lsearch = 'none', 'lsh','lsb' or 'elsr' or 'elsc'
%         - it is the current ALS iteration
%         - data_type = 'real' or 'complex' 
%------------------------------------------------------------------------------------------------------
% OUTPUTS: Interpolated loading matrices ready to be used in the next ALS update
%------------------------------------------------------------------------------------------------------
% Line Search at the n th ALS iteration consists of interpolation of the loading 
% matrices A1, A2, A3, A4 from their estimates at (n-1)th iteration (given by A1p,A2p,A3p,A4p) 
% and at (n-2)th iteration (given by A1q,A2q,A3q,A4q).
% The crucial point is to find a good step size in the interpolation directions
% This step size is denoted STEP here.
%-------------------------------------------------------------------------------------------------------

it_start=2;   % Line Search will start after it_start iterations (at least 2)

% First compute the interpolation directions
GA1=A1p-A1q;
GA2=A2p-A2q;
GA3=A3p-A3q;
GA4=A4p-A4q;

%-----------------------------------------------------
% 'none'            
%-----------------------------------------------------
if   strcmp(lsearch,'none')==1
     A1=A1p;A2=A2p;A3=A3p;A4=A4p;     % Matrices are not interpolated

%-----------------------------------------------------
% 'lsh'            
%-----------------------------------------------------
elseif strcmp(lsearch,'lsh')==1     
     if it<it_start
         A1=A1p;A2=A2p;A3=A3p;A4=A4p;     % Matrices are not interpolated
     else
         % Interpolatation
         STEP=1.25;
         A1=A1q+STEP*GA1;
         A2=A2q+STEP*GA2;
         A3=A3q+STEP*GA3;
         A4=A4q+STEP*GA4;
         % Compute phi after Interpolation
         phiint=norm(X4-kr(A1,A2,A3)*A4.','fro');
         % Accept or Reject Interpolation                   
         if phiint > phi     % reject interpolation
             A1=A1p;A2=A2p;A3=A3p;A4=A4p;
         end
     end
%-----------------------------------------------------
% 'lsb'
%-----------------------------------------------------
elseif strcmp(lsearch,'lsb')==1     
     if it<it_start
         A1=A1p;A2=A2p;A3=A3p;A4=A4p;     % Matrices are not interpolated
     else
         % Interpolatation
         STEP=it^(1/3);
         A1=A1q+STEP*GA1;
         A2=A2q+STEP*GA2;
         A3=A3q+STEP*GA3;
         A4=A4q+STEP*GA4;
         % Compute phi after Interpolation
         phiint=norm(X4-kr(A1,A2,A3)*A4.','fro');
         % Accept or Reject Interpolation                   
         if phiint > phi     % reject interpolation
             A1=A1p;A2=A2p;A3=A3p;A4=A4p;
         end
     end  

%-----------------------------------------------------
% elsr and elsc
%-----------------------------------------------------
elseif (strcmp(lsearch,'elsr')==1) || (strcmp(lsearch,'elsc')==1)     
     if it<it_start
         A1=A1p;A2=A2p;A3=A3p;A4=A4p;     % Matrices are not interpolated
     else
         % Interpolatation
         STEP=ELS_4way(A1q,A2q,A3q,A4q,GA1,GA2,GA3,GA4,X4,lsearch);
         A1=A1q+STEP*GA1;
         A2=A2q+STEP*GA2;
         A3=A3q+STEP*GA3;
         A4=A4q+STEP*GA4;
     end      
end
end

%********************************************************************************   
function [STEP]=ELS_4way(A1,A2,A3,A4,GA1,GA2,GA3,GA4,X4,lsearch,Tol,MaxIt_in)
% Enhanced Line Search for Cp4. 
%-----------------------------------------------------------------------------
% [7]  M. Rajih, P. Comon and R.A. Harshman, "Enhanced Line Search: A Novel Method to Accelerate 
%     PARAFAC", SIAM Journal Matrix Anal. and Appl. (SIMAX), Volume 30 , Issue 3, pp. 1128-1147, Sept 2008 
% [8]: D. Nion and L. De Lathauwer, "An Enhanced Line Search Scheme for Complex-Valued Tensor Decompositions. 
%      Application in DS-CDMA", Signal Processing, vol. 88, issue 3, pp. 749-755, March 2008.
%-------------------------------------------------------------------------------------------------------
if exist('Tol')~=1 || isempty(Tol)    % Tolerance for the inner loop of ELS with complex step
    Tol=1e-4;
end 

if exist('MaxIt_in')~=1 || isempty(MaxIt_in)    % max nb of iterations for the inner loop of ELS with complex step
    MaxIt_in=10;
end

%---- 4 G terms
Mat4 = kr(GA1,GA2,GA3)*GA4.';

%----- 3 G terms
Mat3 = kr(GA1,GA2,GA3)*A4.'...
   + ( kr(GA1,GA2,A3) + kr(GA1,A2,GA3) + kr(A1,GA2,GA3) )*GA4.';

%----- 2 G terms
Mat2 = ( kr(GA1,GA2,A3) + kr(GA1,A2,GA3) + kr(A1,GA2,GA3) )*A4.' ...
   + ( kr(A1,A2,GA3) + kr(A1,GA2,A3) + kr(GA1,A2,A3) )*GA4.';

%------ 1 G term
Mat1 = ( kr(GA1,A2,A3) + kr(A1,GA2,A3) + kr(A1,A2,GA3) )*A4.' ...
   + kr(A1,A2,A3) *GA4.';

%----- 0 G term
Mat0 = kr(A1,A2,A3)*A4.'-X4;

% Build matrix of polynomial coefficients
M=[Mat4(:) Mat3(:) Mat2(:) Mat1(:) Mat0(:)];
H=M'*M;

% Given the coefficients in H, find the optimal STEP
%--------------------------REAL CASE----------------------------------------     
if strcmp(lsearch,'elsr')==1  % choose the best real STEP 

    % The coefficients of the 8th order polynomial correspond to the sum of the different antidiagonal terms of H
    H2=fliplr(H);  % anti-diagonal is now diagonal
    k=0;
    for n=4:-1:-4
        k=k+1;
        pol(1,k)=sum(diag(H2,n));
    end
    pol_der=polyder(pol);
    sqrts=roots(pol_der);
    sqrts=sqrts(imag(sqrts)==0); % real roots
    sqrts=[sqrts;1];   
    extremum=polyval(pol,sqrts);
    STEP=sqrts(find(extremum==min(extremum),1));
    STEP=STEP(1);               

%--------------------------COMPLEX CASE----------------------------------------                             
elseif strcmp(lsearch,'elsc')==1  

    % Init
    STEP=1;
    a=1;  % real part of step
    b=0;  % imag part of step
    u=[STEP^4 STEP^3 STEP^2 STEP 1].';
    norm_new=abs(u'*H*u);
    norm_diff=norm_new;
    it_in=0;

    % use the following code with symbolic variable to get quickly the expressions of (a+i*b)^N
    % syms a;syms b;L=a+j*b;expand(L^4)
    % or simply use Pascal triangle formulas

    % Alternate update of real part a and imag part b of STEP 
    while (norm_diff > Tol) && (it_in < MaxIt_in)
           it_in=it_in+1;
           norm_old=norm_new;
           %---------------------
           % Update imag part b  
           %---------------------
           Da=zeros(5,5);
           Da(1,:)=[1 -4*j*a -6*a^2 4*j*a^3 a^4];  % from (a+jb)^4
           Da(2,2:end)=[-j -3*a 3*i*a^2 a^3];      % from (a+jb)^3
           Da(3,3:end)=[-1 2*i*a a^2];             % from (a+jb)^2
           Da(4,4:end)=[j a];                      % from (a+jb)
           Da(5,5)=1;
           % Given a fixed, compute the coefficients of the 5th order polynomial of variable b
           Ha=Da'*H*Da;
           Ha2=fliplr(Ha);  % sum of anti-diagonal of Ha are the polynomial coefficients
           k=0;
           for n=4:-1:-4
               k=k+1;
               pol(1,k)=sum(diag(Ha2,n));
           end
           pol_der=polyder(pol);
           sqrts=roots(pol_der);
           sqrts=sqrts(imag(sqrts)==0);
           sqrts=[sqrts;1];
           extremum=polyval(pol,sqrts);
           b=sqrts(find(extremum==min(extremum),1));
           b=b(1);

           %---------------------
           % Update real part a  
           %---------------------
           Db=zeros(5,5);
           Db(1,:)=[1 4*j*b -6*b^2 -4*j*b^3 b^4];  % from (a+jb)^4
           Db(2,2:end)=[1 3*j*b -3*b^2 -j*b^3];      % from (a+jb)^3
           Db(3,3:end)=[1 2*i*b -b^2];             % from (a+jb)^2
           Db(4,4:end)=[1 j*b];                      % from (a+jb)
           Db(5,5)=1;
           % Given b fixed, compute the coefficients of the 5th order polynomial of variable a
           Hb=Db'*H*Db;
           Hb2=fliplr(Hb);  % sum of anti-diagonal of Ha are the polynomial coefficients
           k=0;
           for n=4:-1:-4
               k=k+1;
               pol(1,k)=sum(diag(Hb2,n));
           end
           pol_der=polyder(pol);
           sqrts=roots(pol_der);
           sqrts=sqrts(imag(sqrts)==0);
           sqrts=[sqrts;1];
           extremum=polyval(pol,sqrts);
           a=sqrts(find(extremum==min(extremum),1));
           a=a(1);

           %------------
           % update STEP
           %------------
           STEP=a+j*b;
           u=[STEP^4 STEP^3 STEP^2 STEP 1].';
           norm_new=abs(u'*H*u);
           norm_diff=abs(norm_new-norm_old);
    end
end
end
          