 function [rhs2] =  GdotbG(uhat,b,p)
% Approximation to \phi_z (i.e. Dirichelt-to-Neumann operator)
w = b.*ifft2(p.KxiKy.*uhat);
A = fft2(w);
As = conj(A(p.shift1,p.shift2));
rhs2 = p.KxmiKy.*A/2 + p.KxiKy.*As/2;
end