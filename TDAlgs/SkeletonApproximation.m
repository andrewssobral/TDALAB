function [C, G, R] = SkeletonApproximation(A, r)
%
% This function provides a numerical implementation of the Skeleton
% Approximation Algorithm described in [1]. 
%
% Input: 
%    - A: m x n matrix.
%    - r: number of columns and rows to be selected from A.
%
% Ouput:
%    - C: m x r matrix containing r columns from A.
%    - R: r x n matrix containing r rows from A.
%    - G: r x r matrix.
%
% Outcome:
%    The rank-r matrix T = CGR is ''close'' to the input matrix A, i.e. the
%    approximation error ||A - T||_2 is ''close'' to ||A - Ar||_2, where Ar
%    is the ''best'' rank-r matrix computed by the Singular Value
%    Decomposition of A. 
%
% Notes:
%    The authors of [1] give Theorem 3.1 and describe a constructive proof
%    for their Theorem. The implementation that we describe here is a
%    special case of Theorem 3.1; in particular we set F = A - Ar, and the
%    matrices \hat{U} and \hat{V} of the Theorem are constructed by running 
%    the pivoted QR decomposition method of [2]. This column/row selection
%    step can be replaced by any Rank Revealing QR algorithm with similar
%    properties (see, for example, Table 2 in [3]). 
%
% References:
% 
% [1] S. A. Goreinov, E. E. Tyrtyshnikov, and N. L. Zamarashkin. 
%     A theory of pseudoskeleton approximations. 
%     Linear Algebra and its applications, Volume 261, Issues 1-3, August
%     1997, Pages 1-21
%
% [2] Gene Golub.
%     Numerical methods for solving linear least squares problems. 
%     Numerische Mathematik, 7:206–216, 1965.
%
% [3] C. Boutsidis, M. W. Mahoney, and P. Drineas.
%     An improved approximation algorithm for the column subset selection 
%     problem.
%     Proceedings of the Twentieth Annual ACM-SIAM Symposium on Discrete 
%     Algorithms, SODA 2009, New York, NY, USA, January 4-6, 2009. 
%
% -------------------------------------------------------------------------
% Author: Christos Boutsidis, August 2009. 
% E-mail: christos.boutsidis@gmail.com
% -------------------------------------------------------------------------

% The first step of the algorithm (see eqn. 3.2 in [1]) requires the
% computation of the Singular Value Decomposition of the matrix A - F. 
% Since F = A - Ar, it is A - F = Ar, thus we only need to compute the
% top-r singular values and vectors of the matrix A.
[U, S, V] = svds(A,r);

% Construct the matrix C containing r columns from A.
[qV, rV, pV] = qr(V', 0);
C = A(:, pV(1:r));

% Construct the matrix R containing r rows from A.
[qU, rU, pU] = qr(U', 0);
R = A(pU(1:r),:); 

% Construct the matrix G. We set G as the solution of the following problem
%                       min_G || A - CGR ||_F. 
% This choice of G is in general different from the one in [1]. We chose
% this G since we find its construction simpler than the one in [1] and we 
% believe that - in practice - is as ''good'' as the one in [1].
G = pinv(C, .05) * A * pinv(R, .05); 