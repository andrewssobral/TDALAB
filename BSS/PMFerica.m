function [y w]=PMFerica(x,opts)
% call Equivariant Robust ICA
% 
defopts=struct('NumOfComp',[],'MaxIter',1000,'Tol',1e-3,'whitening',false);
if ~exist('opts','var')
    opts=struct;
end
[nsou,MaxIter,stop,whitening]=scanparam(defopts,opts);
M=size(x,2);
if isempty(nsou)
    nsou=M;
end

[z temp1 temp2]=svds(x,nsou,'L');
z=z';
    
y=erica(z,nsou,whitening,MaxIter,stop);
y=y';
w=y\x;


function [y] = erica(x,nsou,whitening,MaxIter,stop)
% ERICA: Equivariant Robust Indepedent Component Analysis algorithm
% (Asymptotically equivariant for sufficient number of samples
%  in presence of Gaussian noise).
%
%  (c) Sergio Cruces, Luis Castedo, Andrzej Cichocki.
% 	http://viento.us.es/~sergio     E-mail: sergio@us.es
%   Version: 2.5                    Last update:   27/12/2002.
% 					                (Version: 2.0, 28/07/2001)
%
%======================================================
%
% PURPOSE:    To separate the signals from a mixture
%			  of 'n' sensors and 'nsou' sources (with non-zero 
%			  kurtosis) in presence of Gaussian noise. 
%			  The separation solution is shown to be 
%	          a saddle point of a given functional for any
%			  source density whose kurtosis is nonzero.
%			  The algorithm is a quasi-Newton iteration 
%			  that will converge to this saddle point 
%			  with local isotropic convergence regarless
%			  of the sources' densities. The use of 
%			  prewhitening is not necesary for the algorithm
%			  to converge.
%
% MANDATORY INPUT ARGUMENTS
%   x  Matrix with the data observations  
%		 with dimension (Number of sources) x (Number of Samples).
%      Each row corresponds to a different sensor.
%
% OPTIONAL INPUT ARGUMENTS
%   nsou      Number of independent componentes to 
%				  extract simultaneously (default e=1).
%   MaxIter   Maximum number of iterations.
%   stop      Sensivity parameter to stop the convergence. 
%
%
% OUTPUT ARGUMENTS
%   BW  		   Composite ICA transformation.
%   B				Semi-orthogonal transformation.
%   W  		   Prewhitening transformation if selected.
%   A				Estimated mixing system.
%	y				Estimated sources.
%
%=========================================================
%
% Related Bibliography:
%
% [1] S. Cruces, L. Castedo, A. Cichocki,
%		"Blind Source Extraction in Gaussian Noise", 
%		Neurocomputing, vol 49, pp 87-118, Dec. 2002.
%
% [2] S. Cruces, L. Castedo, A. Cichocki,
%		IEEE International Conference on Acoustics,
%		Speech, and Signal Processing, vol. V, pp. 3152--3155,
%       Istambul, Turkey, June 2000 (ISBN: 0-7803-6296-9)
%
% [3] S. Cruces, "An unified view of BSS algorithms",
%		PhD. Thesis (in Spanish), University of Vigo, 
%		Spain, 1999.
%
%=========================================================

% PARAMETERS 

[n,T]=size(x);  % n Sensors Number.
				% T Observations length.
eta0=.9;  		% Learning step size parameter.

if ~exist('nsou','var'), nsou=n;end;
if ~exist('whitening','var'), whitening=0; end 
if ~exist('MaxIter','var'), MaxIter=1000; end 
if ~exist('stop','var'), stop=1e-3; end;

% KS
I=eye(nsou); delta=I; 
B=I;

reinic=2;

% PREWHITENING
x=x-mean(x,2)*ones(1,T);
Rxx=x*x'/T;
if whitening, 
  [u,d,v]=svd(Rxx+eye(n));
  d=diag(d)-1; 
  n=find(d>1e-14, 1, 'last' ); 
  W=(u(:,1:n)*diag(real(sqrt(1./d(1:n))))*u(:,1:n)'); 
else
  W=eye(n,n);
end 
x=W*x;   % Prewhiten the data.

[n,T]=size(x); % n Sensors Number.

B = eye(size(W));
if n<nsou 
   B = eye(n,nsou)';
   I = eye(nsou)';
   nsou=n; 
end   

it=0; 
convergence=0;
while ~convergence 
   it=it+1;
   
   % Output
   y=B*x; 
   
   % Sample cumulants
   
   y_=conj(y);              % complex conjugate output y*
     					    % Cumulant(y,y,y*,y*)
        C11=y*y_.'/T;
        %C12=(y*(y.^2)')/T;
        C13=(y*(y_.*y.*y_).'-2*(y*y_.')*diag(mean((y.*y_).'))...
       -(y*y.')*diag(mean((y_.*y_).')))/T; 
    				        % Diagonal matrix of cumulant signs
   S3=sign(real(diag(diag(C13))));    
   %S2=sign(real(diag(diag(C12))));    

   eta=eta0/(1+fix(it/50)/1.1);   % annealing step-size
   
   pp=C13*S3;
   [ppw,ppk]=size(pp);
   [iw,ik]=size(I);
   if ppw~=iw || ppk~=ik,
      %fprintf( '\n .\n\n' );
      %break
   end
   %if it<5 S3=I; % TRY SPARSE FIRST
   %end  
   Delta=(C13*S3-I);  
   mu=min([eta/(1+eta*norm(Delta+I,1)),2*eta0/(1+3)]);
   %B=(I-mu*Delta)*B;                  % Batch implementations of ERICA (CII) 
   B=inv(I+mu*Delta)*B;                % Alternative implementation: CEASI 
   %B=inv(I+mu*Delta/2)*(I-mu*Delta/2)*B;% Combination of both the previous ones.
   %B=2*inv((I+C11)*(I+mu*Delta))*(I-mu*Delta)*B;  % With implicit prewhitening
  
  % Alternative implementation: CEASI 
  %  Delta=(C11-I)+1/2*(C13*S3-(C13*S3)');
  %  mux=eta/(1+eta*norm(Delta+I,1));
  %  mu1=mux; mu2=min([mux,eta/max(abs(diag(C13)))]); 
  %  B=(I-mu1*(C11-I)-mu2/2*(C13*S3-(C13*S3)'))*B;  
   
     
    phi(1,it)=(sum(sum(abs(Delta))))/(nsou*max((nsou-1),1));   
	delta=phi(1,it);
  	convergence =(delta<stop) | (it==MaxIter);
   
   if (trace(B*B')>1e6) && (reinic>=0), 
   	B=rand(nsou,n);  
	   if reinic==0, 
      	B=0; convergence=1; 
   	else
      	reinic=reinic-1;
   	end         
	end;   
    
end



