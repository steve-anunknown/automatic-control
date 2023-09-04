function rootlocus(nominator, denominator, positive)
m = length(nominator)  -1;
n = length(denominator)-1;
if (n < m)
    fprintf("System has to be strictly proper.\n")
    return
end

sys_zeroes   = roots(nominator);
sys_poles    = roots(denominator);
fprintf("\nThe root locus has %d locus/loci and %d asymptotes.\n", m, n-m)
centroid = (sum(sys_poles)-sum(zeros))/(n-m);
fprintf("\nThe centroid of the asymptotes is %d.\n", centroid)
fprintf("\nThe angle(s) of departure of the asymptotes is:\n")
if (positive)
    for index=1:(n-m)
        theta = (2*index - 1)*pi;
        theta = rad2deg(theta);
        fprintf("\t%d)\tangle = %d\n", index, theta)
    end
else
    for index=1:(n-m)
        theta = (2*(index-1))*pi;
        theta = rad2deg(theta);
        fprintf("\t%d)\tangle = %d\n", index, theta)
    end
end


reals   = sort([real(sys_zeroes); real(sys_poles)]);
uniques = unique(reals);
occurs  = [groupcounts(reals); 0];
in      = zeros(length(uniques)+1, 1);
in(end) = false + ~positive;
for index=length(occurs):-1:2
    if (mod(occurs(index-1), 2) == 0)
        in(index-1) =  in(index);
    else
        in(index-1) = ~in(index);
    end
end


fprintf("\nThe parts of the real axis that belong to the root locus are:\n")
if (in(1))
    fprintf("(-inf, ")
end
index = 1;
while (index <= length(in))
    while (index <= length(in) && in(index) == true)
        index = index + 1;
    end
    if (index == length(in) + 1)
        fprintf("+inf)")
    elseif (index-1 > 0 && in(index-1) ~= false)
        fprintf("%d)", uniques(index-1))
    end
    while (index <= length(in) && in(index) == false)
        index = index + 1;
    end
    if (index <= length(in) && in(index) == true)
        fprintf("(%d, ", uniques(index-1))
    end
end
fprintf("\n")
fprintf("\nThe departure angles from the poles are:\n")
occurs = groupcounts(sys_poles);
for index=1:length(occurs)
    angles_pole_zeros = sum(angle(sys_poles(index)-sys_zeroes));
    angles_pole_poles = sum(angle(sys_poles(index)-sys_poles));
    theta = (1/occurs(index))*(-pi*positive+angles_pole_zeros+angles_pole_poles);
    theta = wrapTo360(rad2deg(theta));
    fprintf("\tpole = %d + i*(%d) => angle = %d\n", real(sys_poles(index)), imag(sys_poles(index)), theta)
end
fprintf("\nThe arrival angles at zeroes are:\n")
occurs = groupcounts(sys_zeroes);
for index=1:length(occurs)
    angles_zero_poles  = sum(angle(sys_zeroes(index)-sys_poles));
    angles_zero_zeroes = sum(angle(sys_zeroes(index)-sys_zeroes));
    theta = (1/occurs(index))*(pi*positive+angles_zero_poles-angles_zero_zeroes);
    theta = wrapTo360(rad2deg(theta));
    fprintf("\tzero = %d + i*(%d) => angle = %d\n", real(sys_zeroes(index)), imag(sys_zeroes(index)), theta)
end

fprintf("\nThe possible insertion points are:\n")
[a, ~] = polyder(denominator, nominator);
possible = roots(a);
for index=1:length(possible)
    fprintf("\t%d)\t possible point = %d", index, possible(index))
    if (~isreal(possible(index)))
        fprintf("\trejected because it is complex\n")
    else
        fprintf("\n")
    end
end
fprintf("Check whether they belong in the root locus or not in order to accept or reject them.\n")
end