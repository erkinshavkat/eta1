function R_hat=Rhat_mild(eta0_hat,eta0p_hat,t,p)
    % x , y : position to evaluate R
    % t : current time value
    % xp, yp : position of past impacts
    % H_data, dH_data : (Nk, Nimpacts), H and dH associated with each past impact evaluated at time t
    % p : parameter structure


    [eta0p_hat,eta0pp_hat] = eta0k_rhs(eta0_hat,eta0p_hat,t,p);

    q1= 2*p.nu0.*p.K2_deriv.*eta0_hat - eta0p_hat;
    q2 = 2*p.nu0.*p.K2_deriv.*eta0p_hat - eta0pp_hat;
    
    hatbq1= fft2(p.d1.*real(ifft2(q1)));
    hatbq2= fft2(p.d1.*real(ifft2(q2)));
    R_hat= 1/p.d0.* ...
        ( 2*p.nu0.*p.K2_deriv .* hatbq1-hatbq2);
    % R_func= real(ifft2(R_hat));
end