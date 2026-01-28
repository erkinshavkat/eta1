function integral = eta0_num(x, y, xp,yp,H,p)
    %% Formerly B4, now 3.10, computes eta at (x,y)
    % x, y (Nx x Nx), position to valuate eta
    % xdata, y_data (1 x 1), position of past impacts
    % H_data (Nk, 1), H associated with each past impact
    r = sqrt((x - xp).^2 + (y - yp).^2 );
    oldsize=size(r);
    r= reshape(r,[1 oldsize]);
    integral =  p.b4_prefactor*p.dk * sum(p.K3_vec .* H.*  besselj(0, p.K_vec .* r));
    integral=reshape(integral,oldsize);
end