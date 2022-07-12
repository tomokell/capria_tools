function [U, S, V]  =   lsvd(A, r)
%
%   [U, S, V]   =   lsvd(A, r)
%
%   U, V are the left and right singular vectors of A
%   S contains the singular values
%
%   A is a real or complex input matrix, where size(A,1) > size(A,2)
%   r is the desired output rank of U*S*V' where 1 <= r <= size(A,2)
%
%   Mark Chiew
%   September 2012
%   Updated May 2014
%
%   Greatly speeds up computation of SVD
%   Approximate SVD calcluation to speed up SVDs for large matrices
%   Assumes A has more rows than columns
%   Computes SVD by first computing the eigenvalues and eigenvectors of the
%   short axis product A'*A
%   This computes r squared singular values (i.e. S) and r right singular
%   vectors of A (i.e. V) respectively.
%   The r left singular vectors of A (i.e. U) are then found by a least squares
%   approximation.
%   The resultant product U*S*V' is a rank r approximation of A, with good
%   fidelty compared to the builtin SVD or SVDS.
%   ============================================================================
    if nargin == 1
        r   =   min(size(A));
    end

    [V, D]  =   eig(A'*A);

    S   =   D.^0.5;
    S   =   S(end:-1:end-r+1,end:-1:end-r+1);
    V   =   V(:, end:-1:end-r+1); 
    U   =   A*(V*diag(1./diag(S)));
