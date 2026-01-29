 function [rhs2] =  GdotbG(uhat,b,p)
% Compute \nabla \cdot( b *\nabla u) in F-space
w = b.*ifft2(p.KxiKy.*uhat);
A = fft2(w);
As = conj(A(p.shift1,p.shift2));
rhs2 = p.KxmiKy.*A/2 + p.KxiKy.*As/2;
end