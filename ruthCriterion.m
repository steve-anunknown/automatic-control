function ruthCriterion(flipped)
    fprintf("\n")
    degree = length(flipped);
    if (mod(degree, 2) == 0)
        cols = degree/2;
    else
        cols = ceil(degree/2);
    end
    A = zeros(degree + 1, cols);
    [rows, cols] = size(A);
    A(1, 1:length(flipped(1:2:end))) = flipped(1:2:end);
    A(2, 1:length(flipped(2:2:end))) = flipped(2:2:end);
    fprintf("Routh Table:\n")
    fprintf("s^%d: %s\n", degree-1, join(string(A(1, :)), ', '))
    fprintf("s^%d: %s\n", degree-2, join(string(A(2, :)), ', '))
    flag = false;
    for i = 3:rows-1
        for j = 1:cols-1
            A(i, j) = calc_elem(A(i-2:i-1, :), j);
        end
        if (all(A(i, :) == 0))
            fprintf("s^%d: %s Row of zeroes detected.\n\n", degree-i, join(string(A(i, :)), ', '))
            fprintf("Check roots of helper polynomial.\n\t")
            helper = alternate(A(i-1, :), 0);
            printpoly(helper)
            polyroots = roots(helper);
            fprintf("Roots are %s", join(string(polyroots), ', '))
            fprintf(".\n\n")
            if (not(any(real(polyroots) == 0)))
                fprintf("They have non zero real part, therefore system is unstable.\n")
                return
            else
                fprintf("They are purely imaginary, therefore the test goes on.\n\n")
                if (even(i))
                    fprintf("The line is even, therefore system is unstable (BIBO and Lyapunov).\n")
                else
                    fprintf("The line is odd, therefore use the coefficients \nof the derivative of the helper polynomial in the \nfollowing line.\n\n")
                    b       = alternate(A(i-1, :), 0);
                    bder    = polyder(b);
                    A(i, :) = [bder, zeros(cols-length(bder))];
                    fprintf("Derivative of helper polynomial: ")
                    printpoly(A(i, :))
                    fprintf("\n")
                    flag    = true;
                end
            end
        end
        if (A(i, 1) == 0)
            A(i, 1) = eps;
        end
        fprintf("s^%d: %s\n", degree-i, join(string(A(i, :)), ', '))
    end
    fprintf("\n")
    sign_changes = count_sign_changes(A(:, 1));
    if (sign_changes > 0)
        fprintf("System is unstable with %d poles on the right half plane.\n", sign_changes)
        return
    end
    if (not(flag))
        fprintf("No sign changes in the first column,\ntherefore system is stable.\n")
    else
        polyroots = roots(b);
        if (length(polyroots) == length(unique(polyroots)))
            fprintf("The roots of the helper polynomial are all distinct,\ntherefore system is marginally stable.\n")
        else
            fprintf("System is unstable.\n")
        end
    end
end

function printpoly(coeffs)
    flipped = flip(coeffs);
    for i = length(flipped):-1:1
        if (flipped(i) ~= 0)
            if (i == length(flipped))
                fprintf("%ds^%d", flipped(i), i-1)
            else
                if (flipped(i) > 0)
                    fprintf(" + %ds^%d", flipped(i), i-1)
                else
                    fprintf(" - %ds^%d", -flipped(i), i-1)
                end
            end
        end
    end
    fprintf("\n")
end

function count=count_sign_changes(col)
    count = 0;
    for i = 2:length(col)
        if (col(i) * col(i-1) < 0)
            count = count + 1;
        end
    end
end

function a=calc_elem(table, i)
    mat = [table(1, 1), table(1, i+1);
           table(2, 1), table(2, i+1)];
    a   = -det(mat)/table(2, 1);
end

function bool=even(num)
    if (mod(num, 2) == 0)
        bool = true;
    else
        bool = false;
    end
end

function arr=alternate(arr1, const)
    arr          = zeros(1, 2*length(arr1));
    arr(1:2:end) = arr1;
    arr(2:2:end) = const;
    pointer = length(arr);
    while (arr(pointer) == 0)
        pointer = pointer - 1;
    end
    arr = arr(1:pointer);
end