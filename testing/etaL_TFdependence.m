%% example comparing etaL computed with exact HL, exact HLinf, and the approximate etaL formula
clear; close all; clc
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=128;
Nk=Nx; Lk=12*pi;
mem=0; theta=0;
eta0=zeros(Nx);

p80 = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,80);
p800 = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,800);

nimpacts=100;
xdim=linspace(-0.02,0.02,Nx);
x80=xdim./p80.xF;
x800=xdim./p800.xF;

eta80_ax=plot(xdim,zeros(Nx,1),"LineWidth",2,"DisplayName","80Hz"); hold on
eta800_ax=plot(xdim,zeros(Nx,1),"LineWidth",2,"DisplayName","800Hz"); hold on

% ylim([-0.002 0.002])
legend()
%preallocate HL

impact_times=[];
for tdim=0:0.0001:10
    

    eta80_ax.YData=etaL_approx(x80,0,tdim/p80.TF,p80).*p80.xF;
    eta800_ax.YData=etaL_approx(x800,0,tdim/p800.TF,p800).*p800.xF;
    pause(0.01)

end