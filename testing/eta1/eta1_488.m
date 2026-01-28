clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=128;
Nk=Nx; Lk=12*pi;
mem=0.98; theta=1.3;
eta0=zeros(Nx);
omega=80;
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);

epsilon=gather(p.epsilon);
nimpacts=100;


t = p.theta/(4*pi);

%setup H code
H_vec=zeros(Nk,nimpacts);
dH_vec=zeros(Nk,nimpacts);

%setup Faria code
phi = p.phi0; eta = p.eta0; 
xi=0; yi=0; ui=0; vi=0; 
phi_hat = fft2(phi);eta_hat = fft2(eta);


%plotting
eta1_plot_domain=-4:p.hx:4;
faria_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","Faria"); hold on
eta0_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","\eta_0"); hold on
eta1_ax=plot(eta1_plot_domain,zeros(size(eta1_plot_domain)),"LineWidth",2,"DisplayName","\eta_0+\epsilon\eta_1"); hold on
legend()
xlim([-4 4]);
ylim([-0.02,0.02])
yyaxis right
topo=plot(p.x,p.d_grid(end/2,:),"Color",'black',"LineWidth",1,"DisplayName","topo"); hold on
ylim([0.1,0.25])

xps=zeros(nimpacts,1);
yps=zeros(nimpacts,1);
R_data=zeros(nimpacts*p.nsteps_impact,Nx,Nx);
H_data=zeros(Nk,nimpacts*p.nsteps_impact);
dH_data=zeros(Nk,nimpacts*p.nsteps_impact);
RH_index=1;
for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])
    [~, ~, phi_hat] = drop_impact_matt(xi,yi, ui, vi, phi_hat, eta_hat, p);

    dH_vec(:,n)=1;
    
    for nn=1:p.nsteps_impact 
        t=t+p.dt;
        dH_data(:,RH_index)=1;
        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p); 
        [H_vec, dH_vec] = H_eq_rkstep(H_vec, dH_vec, t, p);
        [H_data,dH_data]=H_eq_rkstep(H_data,dH_data,t,p);
        R_snapshot=R_alt(p.xx,p.yy,t,xps(1:n),yps(1:n),...
           H_vec(:,1:n),dH_vec(:,1:n),p);
        R_data(RH_index,:,:)= R_snapshot;
        RH_index=RH_index+1;

    end
    if n>0
        faria_eta=gather(real(ifft2(eta_hat)));

        eta0=zeros(Nx);
        for impact=1:n
            xp=xps(impact);yp=yps(impact);
            H_impact=H_vec(:,impact);

            eta0=eta0+ eta0_vectorized(p.xx,p.yy,xp,yp,H_impact,p);

        end


        eta1=0*eta1_plot_domain;
        
        for index = 1:length(eta1_plot_domain)
            xpos=eta1_plot_domain(index);ypos=0;
            r=sqrt((xpos-p.xx).^2+(ypos-p.yy).^2);
            eta1_point=0;
            for i=0:RH_index-1
                s=t-i*p.dt;
                Hs=H_data(:,RH_index-i);
                g_arr=g_num(Hs,r,p);
                R_arr=reshape(R_data(RH_index-i,:,:),size(r));
                eta1_point= eta1_point + p.dt/(2*pi).* sum(sum(g_arr.*R_arr.*p.hx.*p.hy));
            end
            eta1(index)=epsilon*eta1_point+eta0_num(xpos,ypos,xps,yps,H_vec,p);
        end
            
        
        eta0=gather(eta0);
        eta0_arr=eta0(Nx/2+1,:);
        faria_ax.YData=faria_eta(Nx/2+1,:);
        eta0_ax.YData=eta0_arr;
        eta1_ax.YData=gather(eta1);
        drawnow;
        pause(0.1)
    end


end