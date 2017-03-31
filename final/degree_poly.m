%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Chinni Chaitanya (EE13B072) and Prafullachandhra (EE16D402)
% Project-1: BCH-encoder-decoder
% EE5160: Error Control Coding
% Name: degree_poly.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = degree_poly(poly_vec)
    l = length(poly_vec);
    while poly_vec(l) == 0
        l = l-1;
    end
    deg = l-1;
end