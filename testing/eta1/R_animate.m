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

eta0=zeros(p.Nx,p.Ny);eta0prime=zeros(p.Nx,p.Ny);
eta0_hat=fft2(eta0);eta0prime_hat=fft2(eta0prime);
%plotting
eta1_plot_domain=-4:p.hx:4;
% faria_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","Faria"); hold on
% eta0_ax=plot(p.x,zeros(Nx,1),"LineWidth",2,"DisplayName","\eta_0"); hold on
% eta1_ax=plot(eta1_plot_domain,zeros(size(eta1_plot_domain)),"LineWidth",2,"DisplayName","\eta_0+\epsilon\eta_1"); hold on


xps=zeros(nimpacts,1);
yps=zeros(nimpacts,1);
R_data=zeros(nimpacts*p.nsteps_impact,Nx,Nx);
H_data=zeros(Nk,nimpacts*p.nsteps_impact);
dH_data=zeros(Nk,nimpacts*p.nsteps_impact);
RH_index=1;

figure;
h_img = imagesc(ones(Nx,Nx));hold on;


for n=1:nimpacts
    
    disp(['Impact number: ' num2str(n)])
    [~, ~, phi_hat] = drop_impact_matt(xi,yi, ui, vi, phi_hat, eta_hat, p);

    [~, ~, eta0prime_hat] = eta0k_impact(xi,yi, ui, vi, eta0prime_hat, p);
    
    for nn=1:p.nsteps_impact 
        t=t+p.dt;
        [phi_hat, eta_hat] = evolve_wave_IF_rkstep(phi_hat, eta_hat, t, p);
        [eta0_hat, eta0prime_hat] = eta0k_rkstep(eta0_hat, eta0prime_hat, t, p);
        R_snapshot = R_num(eta0_hat,eta0prime_hat,t,p);
        R_data(RH_index,:,:)= R_snapshot;
        h_img.CData=gather(R_snapshot);
        colorbar;
        pause(0.1);
        RH_index=RH_index+1;

    end


end