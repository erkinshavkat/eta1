function R_hat=Rhat_num(eta0_hat,eta0p_hat,t,p)
    % x , y : position to evaluate R
    % t : current time value
    % xp, yp : position of past impacts
    % H_data, dH_data : (Nk, Nimpacts), H and dH associated with each past impact evaluated at time t
    % p : parameter structure


    [eta0p_hat,eta0pp_hat] = eta0k_rhs(eta0_hat,eta0p_hat,t,p);

    q1= 2*p.nu0.*eta0_hat - p.Kneg2.*eta0p_hat;
    q2 = 2*p.nu0.*eta0p_hat - p.Kneg2.*eta0pp_hat;
    
    q1_op=GdotbG(q1,p.d1,p);
    q2_op=GdotbG(q2,p.d1,p);
    R_hat= 1/p.d0.* ...
        ( 2*p.nu0.*p.K2_deriv .* q1_op-q2_op);
    % R_func= real(ifft2(R_hat));
end