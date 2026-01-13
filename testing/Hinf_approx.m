clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974;
% H=0.005;
gamma=5.4953;
H=0.001;

Lx=16;
Nx=256;
eta0=zeros(Nx,Nx);
Nk=Nx; Lk=10*pi;
mem=0;
theta=0;
nimpacts=100;
eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk,theta,mem);

alpha = p.K_vec.* sqrt(p.d0 *(p.K_vec.^2*p.Bo + p.G));
alpha_quad = sqrt(p.d0*p.Bo).*p.K_vec.^2;
alpha_quadC = alpha_quad+sqrt(p.d0/p.Bo) *p.G /2;

% loglog(p.K_vec,1./((cosh(2*p.nu0*p.K2_vec)-cos(alpha))),'LineWidth',2,'DisplayName','Full w/ full \alpha');hold on
% loglog(p.K_vec,exp(-2*p.nu0*p.K2_vec)./(2*p.d0*p.G.*p.K2_vec.^2),'LineWidth',2,'DisplayName','Approx exp/k^2');hold on
% % ylim([-10,10])
% legend()
% return
t = p.theta/(4*pi);
t0= p.theta/(4*pi);
t_vec=zeros(nimpacts,1);

[Hinftau0,Hinftau0_approx]=Hinftau(0,p);
plot(p.K_vec/(2*pi),Hinftau0.*p.K3_vec,"LineWidth",2,"DisplayName",'k^3H_L^\infty analytic'); hold on;
plot(p.K_vec/(2*pi),Hinftau0_approx.*p.K3_vec,"LineWidth",2,"DisplayName",'k^3H_L^\infty 451'); hold on;
legend()
saveas(gcf,'vis/line451_800hz','png')
return
HLinf_full_axis=plot(p.x,0*p.x,"LineWidth",2,"DisplayName",'Full H_L^\infty'); hold on;
HLinf_approx_axis=plot(p.x,0*p.x,"LineWidth",2,"DisplayName",'Approx H_L^\infty'); hold on;
% HLinf_full_axis=plot(p.K_vec/(2*pi),0*p.K_vec,"LineWidth",2,"DisplayName",'Full H_L^\infty'); hold on;
% HLinf_approx_axis=plot(p.K_vec/(2*pi),0*p.K_vec,"LineWidth",2,"DisplayName",'Approx H_L^\infty'); hold on;
eta_data=zeros(Nx,p.nsteps_impact);


etaL_Data=zeros(Nx,p.nsteps_impact);
for i=1:nimpacts
    disp(['Impact ' num2str(i)])
    dH_vec(:,i)=1;
    for ii=1:p.nsteps_impact

        etaL_fullHLinf=zeros(Nx,1);
        etaL_approxHLinf=zeros(Nx,1);

        [HLsum,HLapprox]=Hinftau(t-floor(t),p);
        for j=1:Nx
            etaL_fullHLinf(j)=b4_eta(p.x(j),0,0,0,HLsum,p);
            etaL_approxHLinf(j)=b4_eta(p.x(j),0,0,0,HLapprox,p);
        end
        HLinf_full_axis.YData=etaL_fullHLinf;
        HLinf_approx_axis.YData=etaL_approxHLinf;
        ylim([-0.2,0.2])

        % HLinf_full_axis.YData=HLsum.*p.K3_vec;
        % HLinf_approx_axis.YData=HLapprox.*p.K3_vec;
        % ylim([-50,50])
        % xlim([0,5])
        legend()
        pause(0.1)
        t= t+p.dt;

    end
end