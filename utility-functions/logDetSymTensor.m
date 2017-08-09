function ld = logDetSymTensor(A)
% FORMAT ld = logDet(A)
% A  - A field of sparse symetric matrices
% ld - Logarithm of determinant of A
%
% Log-determinant of a sparse squre matrix which has the form of a tensor
% field.
% ld = sum(log(pointwiseDet(A)))

    % Original dimensions of A
    dim = size(A);
    
    % "Reshaper" of the symmetric matrices
    ind = symIndices(dim(end));
    n = size(ind, 1);
    
    % Insure A has 4 dimensions (and thus H has 5 dimensions)
    A = shiftdim(A, length(dim) - 4);
    sdim = size(A);
    
    % Allocate full tensor field
    H = zeros([sdim(1:3) n n], 'single');
    
    % Fill it
    for i1=1:n
        for i2=1:n
            H(:,:,:,i1,i2) = A(:,:,:,ind(i1,i2));
        end
    end
    
    % Compute pointwise determinant
    H = spm_diffeo('det', H);
    
    % Give back original dimensions
    H = shiftdim(H, 4 - length(dim));
    
    % Compute full log determinant
    ld = sumall(log(H(H>0)));

end