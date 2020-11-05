function D_w = Whitening(Wtr,Ntrain,Lr)
    % obtain noise covariance matrix 
    Wtr = Wtr(:,1:Lr*Ntrain);
    blocks = reshape(Wtr, 64,4,Ntrain);
    blocks_l = num2cell(blocks,[1,2]);
    blocks_l = reshape(blocks_l,1,Ntrain);
    blocks_c = cellfun(@(x) x'*x, blocks_l, 'UniformOutput', false);
    C = blkdiag(blocks_c{:});
    D_w = chol(C);
end
