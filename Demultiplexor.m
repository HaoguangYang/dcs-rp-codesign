% Demultiplxing the input matrix K into nBlocks with equal sizes
% (rounding towards the nearest row).
%
% input:
% K: the input matrix to be demultiplexed
% nBlock: number of sub-matrices to divide it into.
%
% output:
% RowBlock: block-ified sub-matrices by row

function RowBlock = Demultiplexor(K, nBlock)
	[KRow, ~] = size(K);
	unit_size = KRow/nBlock;
    RowBlock = cell(nBlock, 1);
    for i = 1:nBlock
		UBound = round(i*unit_size);
		LBound = round((i-1)*unit_size);
		RowBlock{i, 1} = K(LBound+1:UBound,:);
    end
end
