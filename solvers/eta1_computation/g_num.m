function integral=g_num(H_vec,r,p)
    %H_vec (Nk,1) values of H in k-space at given time
    %r 1x1 stands for |\bx - \by|
    oldsize = size(r);
    r= reshape(r,[1 oldsize]);
    integral =  p.dk * sum(p.K_vec .* H_vec.*besselj(0, p.K_vec .* r));
    integral=reshape(integral,oldsize);
end