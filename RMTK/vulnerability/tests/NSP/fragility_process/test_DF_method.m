clear all
clc
close all
addpath('./functions/');
g=9.81;
Tc=0.5;
Td=1.8;
mc = 2.5;
a = 0.00;
ac = 2.36;
r = 0.44;
mf = 6.87;
Gamma =	1.292;				
T	= 1.612;				
Tav	= 1.612;				
Sa_ratios =	1;				
SPO	=  [0.0422728605913 0.113278436093 0.123312633118 0.290448402 2133.26 2133.26 939.913228223]; %d and F
dcroof	= [0.0781	0.1323	0.1869	0.2448	0.2904];
bUth	= [0	0	0.1	0.2	0.3];
dry = SPO(1);
mu = dcroof/dry;
SaT = [];
b = [];
for i = 1:length(dcroof)
    drlim = dcroof(i);
    buthd = bUth(i);
    [Sa50,bTSa,bRSa]=DFfragility(dry,drlim,buthd,mc,r,T,Gamma,Tc,Td,g);
    SaT = [SaT Sa50];
    b = [b bTSa];
end

SaT
b