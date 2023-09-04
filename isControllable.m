function bool = isControllable(A, B)
    C = controlMatrix(A,B);
    bool = (rank(C) == length(A));
end