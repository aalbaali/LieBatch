function [ blkdiag_mat, mat_cell] = blkdiag3d( mat_3d)
    % Converts 3D array into single matrix by putting all the matrices in
    % block-diagonal matrix. E.g., converts [ 2 x 2 x 5] array to matrix of size
    % [ 5 * 2 x 5 * 2].
    % Optionally returns the cell array matrix
    
    % Number of matrices
    K = size( mat_3d, 3);
    mat_cell = mat2cell( mat_3d, size( mat_3d, 1), size( mat_3d, 2), ones(1, K));
    blkdiag_mat = blkdiag( mat_cell{ :});
end