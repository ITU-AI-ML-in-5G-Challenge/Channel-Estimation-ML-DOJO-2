function D_w = Whitening(Wtr,Ntrain,Lr)
    % obtain noise covariance matrix 
    blocks = reshape(Wtr, 64,4,100);
    blocks_l = num2cell(blocks,[1,2]);
    blocks_l = reshape(blocks_l,1,100);
    blocks_c = cellfun(@(x) x'*x, blocks_l, 'UniformOutput', false);
    C = blkdiag(blocks_c{:});
    D_w = chol(C);
end
