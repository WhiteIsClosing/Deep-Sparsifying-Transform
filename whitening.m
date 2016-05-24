function [ Y ] = whitening(Y)
    sigma =Y*Y'/size(Y,2);
    [U,S,V]=svd(sigma);
    epsilon=1E-6;
    Y = U*diag(1./sqrt(diag(S)+epsilon)) * U' * Y;
end

