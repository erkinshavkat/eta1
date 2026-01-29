clear; close all; clc
format long
addpath(genpath(pwd))
% gamma=4.4974; H=0.005;

gamma=5.4953;H=0.001;


Lx=16; Nx=256;
Nk=Nx; Lk=12*pi;
mem=0.98; theta=1.3;
eta0=zeros(Nx);
for omega=[80,160,400,800]

    p = setup_IF_matt(gamma,H,eta0,Nx,Lx,Nk,Lk,theta,mem,omega);
    disp(p.xF)
    disp(p.kC/(2*pi))
end