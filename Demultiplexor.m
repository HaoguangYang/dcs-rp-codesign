% Demultiplxing the input matrix K into nBlocks with equal sizes
% (rounding towards the nearest row).
%
% input:
% K: the input matrix to be demultiplexed
% nBlock: number of sub-matrices to divide it into.
%
% output:
% RowBlock: block-ified sub-matrices by row
% split: number of rows in each sub-matrices

function [RowBlock, split] = Demultiplexor(K, nBlock)
	[KRow, KCol] = size(K);
	unit_size = KRow/nBlock;
	RowBlock = zeros(nBlock, KCol*ceil(unit_size));
	split = zeros(1,nBlock);
    for i = 1:nBlock
		UBound = round(i*unit_size);
		LBound = round((i-1)*unit_size);
		split(i) = UBound-LBound;
		RowBlock(i,:) = reshape(K(LBound+1:UBound,:),1,[]);
    end
end
