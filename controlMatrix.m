function C = controlMatrix(A,B)
    [n, ~] = size(A);
    [~, r] = size(B);
    C = zeros(n, n * r);
    for i=0:n-1
        C(:, (i*r+1):(i*r+r)) = A^i * B;
    end
end