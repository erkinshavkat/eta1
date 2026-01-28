function [eta_hat,etaprime_hat] = eta0k_rkstep(eta_hat,etaprime_hat, t_in, p)
% Evolves wave field in F-space through nsteps of size dt



t = t_in;    


[rhs_eta_1, rhs_etaprime_1] = eta0k_rhs(eta_hat, ...
                                      etaprime_hat, ...
                                      t, p);

[rhs_eta_2, rhs_etaprime_2] = eta0k_rhs(eta_hat      + p.dt/2 * rhs_eta_1, ...
                                          etaprime_hat + p.dt/2 * rhs_etaprime_1, ...
                                          t+p.dt/2, p);

[rhs_eta_3, rhs_etaprime_3] = eta0k_rhs(eta_hat      + p.dt/2 * rhs_eta_2, ...
                                          etaprime_hat + p.dt/2 * rhs_etaprime_2, ...
                                          t+p.dt/2, p);


[rhs_eta_4, rhs_etaprime_4] = eta0k_rhs(eta_hat      + p.dt * rhs_eta_3, ...
                                          etaprime_hat + p.dt * rhs_etaprime_3, ...
                                          t+p.dt, p);


eta_hat      = eta_hat      + p.dt / 6 * (rhs_eta_1 + 2*rhs_eta_2 + 2*rhs_eta_3 + rhs_eta_4);
etaprime_hat = etaprime_hat + p.dt / 6 * (rhs_etaprime_1 + 2*rhs_etaprime_2 + 2*rhs_etaprime_3 + rhs_etaprime_4);


end

     
