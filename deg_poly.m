function [deg,i] = deg_poly(a)
j = 1;
while (j>0 & j<=127)
    if a(j) == 1
        deg = 127 - j;
        i = j;
        j = 0;
    else
        j = j+1;
    end
end

end
        