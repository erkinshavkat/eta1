function integral=g_num(t,s,r,p)
    %H_vec (Nk,1) values of H in k-space at given time
    %r 1x1 stands for |\bx - \by|
    current_t=s;
    H_vec= zeros(length(p.K_vec),1);
    dH_vec= ones(length(p.K_vec),1);
    while current_t<t
        [H_vec, dH_vec] = H_eq_rkstep(H_vec, dH_vec, current_t, p);
        current_t=current_t+p.dt;
    end
    r= reshape(r,[1 size(r)]);
    integral =  p.dk * sum(p.K_vec .* H_vec.*besselj(0, p.K_vec .* r));

end