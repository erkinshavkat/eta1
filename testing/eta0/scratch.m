clear; close all; clc
addpath(genpath(pwd))
gamma=4.4974;
H=0.005;

Lx=16;
Nx=256;
eta0=zeros(Nx,Nx);
Nk=Nx; Lk=6*pi;

nimpacts=10;
mem=1;
theta=0;
t0=theta/(4*pi);
t=theta/(4*pi);
fig = figure('Position', [0, 0, 800, 600]); 
eta0=zeros(Nx);
p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,0,Lk,theta,mem);

alpha= p.K_vec.* sqrt(p.d0 *(p.K_vec.^2*p.Bo + p.G));
plot(p.K_vec,1./(alpha.*(cosh(2*p.nu0*p.K2_vec)-cos(alpha))),'LineWidth',2); hold on
plot(p.K_vec,exp(-2*p.nu0*p.K2_vec)./alpha,'LineWidth',2);hold off
legend('cosh','exp')
return
x_domain=linspace(-5,5,Nx/2);
ax=plot(x_domain,0*x_domain,"LineWidth",2,"DisplayName",'Faria'); hold on
ax2=plot(x_domain,0*x_domain,"LineWidth",2,"DisplayName",'Milewski'); hold on

% v = VideoWriter(sprintf('vis/Faria vs milewski inf limit 800hz strobe.mp4',p.theta/pi),'MPEG-4');
% v.FrameRate = 6;
% open(v);

for t=0.99:100.99 %1:nimpacts*p.nsteps_impact
    impacts=floor(t);
    HL_Faria=zeros(Nk,1);
    for i=0:impacts
        HL_Faria = HL_Faria+HL_varDtN(t-i,p.K_vec,p);
    end
    etaL_Faria=zeros(Nx/2,1);
    for i=1:length(x_domain)
        x = x_domain(i);
        etaL_Faria(i) = eta_varDtN(x,0,0,0,HL_Faria,p);
    end
    ax.YData = etaL_Faria;
    title(['t=' num2str(t),'T_F'])
    ylim([-0.05,0.05])
    xlim([-5,5])
    legend()
    % frame = getframe(gcf);
    % writeVideo(v, frame);
    pause(1/24)

end
close(v);

% HLinf = @(k) sin(p.alpha_full(k))./(2.*p.alpha_full(k).*(cosh(2*k.^2.*p.nu0)-cos(p.alpha_full(k))));
% HLinf_sechapprox = @(k) sin(p.alpha_quadC(k))./(2.*p.alpha_quad(k).*(cosh(2*k.^2.*p.nu0)));

% HLinf_seriesapprox = @(k) sin(p.alpha_quadC(k))./(2.*p.alpha_quad(k)).*(exp(-(2*k.^2.*p.nu0)));

% x_domain = linspace(-5,5,Nx);

% % plot(p.K_vec/(2*pi),HLinf(p.K_vec),"LineWidth",2,"DisplayName",'full H_L^\infty'); hold on
% % plot(p.K_vec/(2*pi),HLinf_sechapprox(p.K_vec),"LineWidth",2,"DisplayName",'k sin(k^2+C)sech(2\nu k^2)'); hold on
% % plot(p.K_vec/(2*pi),HLinf_seriesapprox(p.K_vec),"LineWidth",2,"DisplayName",'k sin(k^2+C)(exp(-2\nu k^2))'); hold on

% plot(x_domain,arrayfun(@(x) b4_eta(x,0,0,0,HLinf(p.K_vec),p),x_domain),"LineWidth",2,"DisplayName",'full H_L^\infty'); hold on
% plot(x_domain,arrayfun(@(x) b4_eta(x,0,0,0,HLinf_sechapprox(p.K_vec),p),x_domain),"LineWidth",2,"DisplayName",'H_L~k sin(k^2+C)sech(2\nu k^2)'); hold on
