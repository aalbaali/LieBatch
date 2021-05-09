function blkdiags = getBlockDiagonals(A, blk_size)
%GETBLOCKDIAGONALS Extract block diagonals from a square matrix.
%
% PARAMETERS
% ----------
% A : [ N x N ] double
%     The matrix from which to extract block diagonals.
% blk_size : double
%     The size of the blocks to extract.
%
% RETURNS
% -------
% blkdiags : [ blk_size x blk_size x length(A) / blk_size ] double
%     The block diagonals, in a 3D array. 
% -------------------------------------------------------------------------
    if isempty(A)
        blkdiags = [];
        return
    end

    if ~mod(blk_size, 1) == 0 || blk_size <= 0
        error("Block size must be positive integer")
    end
    if blk_size > length(A)
        error("Block size too large for this matrix")
    end
    if mod(length(A), blk_size) ~= 0
        error("This block size will not produce an even number of blocks")
    end
    if size(A, 1) ~= size(A, 2)
        error("This function only supports block extraction from square matrices")
    end

    % Extract block diagonals using logical mask
    tf_mask = logical(kron(speye(length(A) / blk_size), ones(blk_size)));
    blkdiags = reshape(full(A(tf_mask)), blk_size, blk_size, []);
end

