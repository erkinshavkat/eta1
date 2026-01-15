function R_func=R_num(eta0,eta0m1,eta0m2,p)

    lap_eta=real(ifft2(p.K2_deriv.*fft2(eta0)));
    lap_etam1=real(ifft2(p.K2_deriv.*fft2(eta0m1)));
    dt_eta= (eta0 - eta0m1)/p.dt;
    dt_etam1= (eta0m1 - eta0m2)/p.dt;

    bu=p.d1.*(2*p.nu0*lap_eta-dt_eta);
    bu_m1=p.d1.*(2*p.nu0*lap_etam1 - dt_etam1);

    lap_bu=real(ifft2(p.K2_deriv.*fft2(bu)));
    dt_bu= (bu - bu_m1)/p.dt;

    R_func=(2*p.nu0*lap_bu-dt_bu)./p.d0;


end