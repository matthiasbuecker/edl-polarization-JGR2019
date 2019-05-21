% ZMM   Computation of magnitude and phase of the complex impedace Z of the
%   one-dimensional Marshall-Madden model.
%
%   Matthias Buecker, May 2019
%
%   Input parameters
%   T:      Temperature [K]
%   C:      Bulk ion concentration [mol/m?]
%   w:      Angular frequency [1/s]
%   L1:     Length of pore 1
%   mp1:    Cation mobility in pore 1 [m^2/V/s]
%   mn1:    Anion mobility in pore 1 [m^2/V/s]
%   L2:     Length of pore 2
%   mp2:    Cation mobility in pore 2 [m^2/V/s]
%   mn2:    Anion mobility in pore 2 [m^2/V/s]

function [absZ,phiZ]=ZMM(T,C,w,L1,mp1,mn1,L2,mp2,mn2)
% Constants
kB=8.617e-5; % Boltzmann's constant kB = 8,617343? 10-5 eV/K
F=96485.3; % Faraday's constant [C/mol]

% Derived parameters
Dp1=mp1*kB*T; % Cation diffusion coefficient in pore 1 [m^2/s]
Dn1=mn1*kB*T; % Anion diffusion coefficient in pore 1 [m^2/s]
Dp2=mp2*kB*T; % Cation diffusion coefficient in pore 2 [m^2/s]
Dn2=mn2*kB*T; % Anion diffusion coefficient in pore 2 [m^2/s]
% Transport numbers
tp1=mp1./(mp1+mn1); % cations pore 1 [-]
tp2=mp2./(mp2+mn2); % cations pore 2 [-]
tn1=1-tp1; % anions pore 1 [-]
tn2=1-tp2; % anions pore 2 [-]
% Time constants tau...
tau1=1./(8*Dp1.*tn1).*L1.^2; % ... in pore 1
tau2=1./(8*Dp2.*tn2).*L2.^2; % ... in pore 2

% Marshall-Madden impedance
% DC impedance rho0
rho0=kB*T/(C*F)*(L1./(Dp1+Dn1)+L2./(Dp2+Dn2)+...
    8*(tn2.*tp1-tn1.*tp2).^2./(L1./tau1+L2./tau2));
% Chargeability eta0
eta0=kB*T/(C*F)*8*(tn2.*tp1-tn1.*tp2).^2./(L1./tau1+L2./tau2)./rho0;
% Impedance Z
Z=rho0.*(1-eta0.*(1-(L1./tau1+L2./tau2)./...
    (L1./tau1.*sqrt(1i*w.*tau1).*coth(sqrt(1i*w.*tau1))+...
    L2./tau2.*sqrt(1i*w.*tau2).*coth(sqrt(1i*w.*tau2)))));

% Computation of output variables
% Impedance magnitude [Ohm-m]
absZ=abs(Z); 
% Phase of impedance [rad]
phiZ=angle(Z); 