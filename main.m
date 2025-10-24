clear all
close all
clc
%% Offline
Model%specify model
ROM%reduced order model
Peak%compute peak-to-peak gain
Filter%define filter and compute peak-to-peak
save('offline')
%% nominal ROM-based MPC
u_max = 100;z_max = 1;%constraints
R = 1e-3;%input cost
x_0r=zeros(n_r,1);delta_0=0;%initial condition
w_max = 0;%no disturbance
%1. solve nominal ROM-based OCP
T_s=0.1;%here: fine resolution to see smooth plot
H_nom=30/T_s;%discrete horizon
MPC_nominal
%nominal has short horizon, since nothing happens later,
%small sampling time only for nicer plots;
%--> comparison computation itme only meaningful wih same T_s,H!
%% solve robust ROM-based MPC
T_s=2;
H = 300/T_s;%discrete horizon
MPC_robust
%% Compare tightening
T_s_validate=1e-1;%set finer time resolution for simulation and validation
Tightening