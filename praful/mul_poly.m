function mul_vec = mul_poly(c, vec, deg, M, PRIM_POLY)
    a_vec = vec(1: deg+1);
    b_vec = gf(zeros(1, c), M, PRIM_POLY);
    c_vec = gf(zeros(1, 84-c-deg-1), M, PRIM_POLY);
    
    mul_vec = [b_vec, a_vec, c_vec];
end