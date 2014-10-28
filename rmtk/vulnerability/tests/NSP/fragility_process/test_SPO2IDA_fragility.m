clc
clear all
close all
%Input data
addpath('./functions/');
Gamma = 1.292;
T = 1.612;
Tav = 1.612;
Sa_ratios = 1;
SPO =[0.0422728605913 0.113278436093 0.123312633118 0.290448402 2133.26 2133.26 939.913228223];
dcroof = [0.0781	0.1323	0.1869	0.2448	0.2904];
noBlg = 0;
w = 1;
bUthd = [0.	0.	0.1	0.2	0.3];
mf = 6.87;
g = 9.81;
MC = 10;
m = dlmread('inputs/idacm.tcl');
r = dlmread('inputs/idacr.tcl');
for j=1:3
    idac(j).m = m(j,:);
    idac(j).r = r(j,:);
end
dry = SPO(1);
Sa = [];
bSa = [];

for i = 1:length(dcroof)
    drlim = dcroof(i);
    buthd = bUthd(i);
    [SaR50,bRSa,SaT50,bTSa]=SPO2IDAfragility(idac, dry,drlim,buthd,mf,T,Gamma,g,MC);
    Sa = [Sa SaT50];
    bSa = [bSa bTSa];
end

Sa
bSa