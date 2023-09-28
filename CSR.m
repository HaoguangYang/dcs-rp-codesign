% Generates the Compressed Sparse Row format of a sparse matrix.
% input:
% A: the input sparse matrix
%
% outputs:
% ptr: pointer to the start of a row
% col: which coloumn the non-zero element is in
% v: vector of non-zero elements in A.

function [ptr, col, v] = CSR(A)
    % find in Transpose for good sequencing to format into CSR.
    [col, row, v] = find(A.');
    ptr = zeros(size(A,1),1);
    j = 1;
    for i = 1:(size(v,1)-1)
        if row(i)~=row(i+1)
            ptr(j) = i;
            j = j+1;
        end
    end
    ptr(end) = size(v,1);
end
