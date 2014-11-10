function [Sa50,bTSa,bRSa]=DFfragility(dry,drlim,buthd,mc,r,T,Gamma,Tc,Td,g)
% Estimate fragility parameters for an MDOF system approximated by an equivalent
% SDOF and using the R-mu-T relationships of Dolsek & Fajfar 2004 paper 
% (Earthquake Engng Struct. Dyn. 2004; 33:1395–1416). Some application details
% and defaults are taken from Dolsek & Fajfar 2005 paper (Earthquake Engng
% Struct. Dyn. 2005; 34:49–66). Dispersion results are not provided and thus
% taken from Ruiz-Garcia & Miranda, assuming elastoplastic backbone. This
% probably underestimates dispersion beyond the capping point, but that is the
% best we can do without going to SPO2IDA. Original D&F paper suggested cov=0.7
% in short period range and 0.4 for medium/long but this is too coarse for our
% needs.
%
% This R-mu-T is suggested by the authors as conservative, since it is
% not based on the median but on the mean mu given R. It is also based on
% Newmark-Hall type spectra. In the following I will try to correct it assuming
% lognormality and the dispersion of Ruiz-Garcia & Miranda (2007). 
%
% Shape of the backbone is elastic-plastic-negative-residual. 
% Parameters
% Dolsek <=> SPO2IDA
%     ms <=> mc (capping point, i.e., end-of-plastic-plateau, ductility)
%     mu <=> mr (end-of-negative-stiffness ductility, ignored by the model)
%     ru <=> rp (residual strength, normalized by the yield strength)
%____________________________________________________________________________
% INPUT
% Cy     : base shear coefficient at yield (presently covered by Gamma and dry)
%          Cy = Say/g = Vy/W, where Say=Sa @ yield, Vy = base shear @ yield
% drlim  : median roof displacement value that defines the fragility. It can result 
%          from any EDP but it should be expressed in terms of the (median)
%          corresponding roof disp (just like we always do in pushovers anyway)
% dry    : roof yield displacement
% buthd  : dispersion (std of log data) characterizing the lognormal
%          distribution of "limit-state roof drift capacity" around drlim
% T      : ESDOF period (sec)
% Gamma  : the participation factor for the roof displacement(>1). Note that for
%          a tall building (or higher-mode influenced one) it is best if you get
%          this value from Say versus droofy estimated from modal response
%          spectrum analysis to include multiple modes. If you use the Gamma of
%          the first one only you will not do well enough (it tends to be
%          lower). See work of Katsanos & Vamva (2014). 
% Tc,Td  : constant accel-constant velocity  and constant velocity-constant
%          displacement corner periods of a Newmark-Hall type spectrum. Default
%          values (roughly taken from Dolsek & Fafjar's (2005) third set of
%          ground motions), are [0.5,1.8].
% mc     : end-of-plastic-plateau capping ductility.
% r      : residual plateau strength divided by yield strength
%
% g      : the value of "g" in units compatible with T and drlim, dy.
%          The default is 9.81m/s2, assuming that dy and drlim are in meters.
% OUTPUT
% Sa50   : the median Sa for the fragility, in units of "g"
% bTSa   : the total dispersion, including capacity and demand dispersions.
% bRSa   : the dispersion due to record-to-record variability (i.e. demand only)
%_____________________________________________________________________________
%
% CREATED 07/Sep/2014 by D.Vamvatsikos
% 
% Example
% T=1; Gamma=1.3; 
% drlim=0.3; dry=0.1; %(in meters)
% buthd=0.2;
% mc=2; r=0.5;
% [Sa50,bTSa]=DFfragility(dry,drlim,buthd,mc,r,T,Gamma);

if nargin<10, g=9.81; end
if nargin<9, Tc=0.5; end
if nargin<8, Td=1.8; end
if nargin<7, error('need at least 7 arguments'); end
if r>1 || r<0.25, error('rp should be within [0.25,1]'); end
if mc<1 || mc>2.5, error('mc should be within [1,2.5]'); end


mlim=drlim/dry;

% From the graphs the authors do not show data for mean ductilities larger 
% than 10 (Fig. 14b). Let's adopt a maximum for the mean ductility of 10.
if mlim>10, error('mlim should be less than 10'); end 

% For a given period, the model is bilinear (piecewise linear with 2 segments)
% in the post-yield range. 
% In other words it provides a trilinear approximation for the mean IDA (mean mu
% given R). The point where this change happens is [mc,Rmc].
% So the mean IDA is defined by 4 points total and 3 segments (incl. elastic):
% Elastic segment: [0,0] -> [1,1]  
%                  Note though the improvement suggested in D&F (2006) that is
%                  not adopted here.
% Capping segment: [1,1] -> [mc,Rmc]
% Post-capping   : [mc,Rmc] -> [mf,Rmf]
%                  where the final point is selected by some arbitrary
%                  Rmf>Rmc, mf=function of Rmf.
% Thus, to get the mean IDA, I only need mf and Rmf (in addition to mc,Rmc).



% Get the (mc50, Rmc) point
Tdstar=Td*sqrt(2-r);
if T<=Tc
   Rmc=0.7*(T/Tc)*(mc-1)+1;
elseif T<=Tdstar
   DT=(T-Tc)/(Tdstar-Tc);
   Rmc=(0.7+0.3*DT)*(mc-1)+1;
else
   Rmc=mc;
end

bthd_mc=1.957*(1/5.876+1./(11.749*(T+0.1))).*(1-exp(-0.739*(Rmc-1)));
mc50=mc/exp(0.5*bthd_mc^2);
mc16=mc50*exp(-bthd_mc);
mc84=mc50*exp(bthd_mc);

% Now set Rmf as twice Rmc and get [mf50,Rmf].

Rmf=2*Rmc;

R0=Rmc;
m0=mc;

if T<=Tc
      c=0.7*sqrt(r)*(T/Tc)^(1/sqrt(r));
elseif T<=Tdstar
      c=0.7*sqrt(r)*(1-DT)+DT;
else
   c=1;
end



% This is the mean mu given R.
mf_mean=(Rmf-R0)/c + m0;
bthd_mf=1.957*(1/5.876+1./(11.749*(T+0.1))).*(1-exp(-0.739*(Rmf-1)));
% For lognormal: Median = mean * exp(0.5*sigma^2)
mf50=mf_mean/exp(0.5*bthd_mf^2);
mf16=mf50*exp(-bthd_mf);
mf84=mf50*exp(bthd_mf);


R=[0,1,Rmc,Rmf];
mu16=[0,1,mc16,mf16];
mu50=[0,1,mc50,mf50];
mu84=[0,1,mc84,mf84];
% plot(mu16,R,mu50,R,mu84,R); grid on


%b=ln(mu50)/lnR, R at 85% of Rlim.



% Now we have fractiles and can invert. Use linear extrapolation for mlim values
% larger than mf. 
R16=interp1(mu84,R,mlim,'linear','extrap');
R50=interp1(mu50,R,mlim,'linear','extrap');
R84=interp1(mu16,R,mlim,'linear','extrap');

% Go to an R value at 85% of the R50 for a biased estimate of "b"
Rlim=0.85*R50;
m50Rlim=interp1(R,mu50,Rlim,'linear','extrap');
b=log(m50Rlim)/log(Rlim);

bRSa=0.5*(log(R84)-log(R16));
bUSa=buthd/b;


CR50 = mlim/R50;
Sa50 = 4*pi^2*drlim./(g*Gamma.*CR50.*T.^2);
bTSa = sqrt(bUSa.^2 + bRSa.^2);

