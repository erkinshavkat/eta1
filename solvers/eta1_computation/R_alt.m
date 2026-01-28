function R_func=R_alt(x,y,t,xp,yp,H_data,dH_data,p)
    % x , y : position to evaluate R
    % t : current time value
    % xp, yp : position of past impacts
    % H_data, dH_data : (Nk, Nimpacts), H and dH associated with each past impact evaluated at time t
    % p : parameter structure

    domsize=size(x);
    u=zeros(domsize);
    Lapu=zeros(domsize);
    dtu=zeros(domsize);

    impacts=size(xp,1);

    [dH_data,d2H_data]=H_rhs(H_data,dH_data, t, p);

    u_kernal = -2*p.nu0*p.K2_vec .* H_data - dH_data;
    Lapu_kernal = -p.K2_vec.* u_kernal;
    dtu_kernal =-2*p.nu0*p.K2_vec .* dH_data - d2H_data;

    for impact=1:impacts
        u=u+ eta0_vectorized(x,y,xp(impact),yp(impact),u_kernal(:,impact),p);
        Lapu=Lapu+ eta0_vectorized(x,y,xp(impact),yp(impact),Lapu_kernal(:,impact),p);
        dtu=dtu+ eta0_vectorized(x,y,xp(impact),yp(impact),dtu_kernal(:,impact),p);

    end
    R_func= (2*p.nu0.*(p.d1.*Lapu+p.lap_d1.*u)-p.d1.*dtu)/p.d0;
end