m_z = zeros(1,49);
m_1 = [m,m_z];
[deg, i] = deg_poly(m_1);
a = m_1;
g = [1,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,0,1];

while ( deg >= 49)
    z1 = zeros(1,i-1);
    z2 = zeros(1,78-i);
    g_1 = [z1,g,z2];
    a = bitxor(a,g_1);
   [deg,i] = deg_poly(a);
end

c = bitxor(m_1,a);
