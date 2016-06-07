function [ Y ,ZCA] = whitening(Y)
    sigma =Y*Y'/size(Y,2);
    [U,S,V]=svd(sigma);
    epsilon=1E-6;
    ZCA=U*diag(1./sqrt(diag(S)+epsilon)) * U';
    Y = ZCA * Y;
end

