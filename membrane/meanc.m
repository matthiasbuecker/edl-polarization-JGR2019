% meanc   Computation of the average cation and anion concentrations for a
%   cylindrical pore cross section and small zeta potentials (f0<25mV)
%   
%   Matthias Buecker, May 2019
% 
%   INPUT PARAMETERS
%   r0:     Pore radius [m]
%   C0:     Ion concentration in the bulk electrolyte [mol/m^3]
%   T:      Temperature [K]
%   f0:     Zeta potential [V]
%   fQ:     Partition coefficient (Leroy et al., 2008) [-]

function [cp_mean,cn_mean]=meanc(r0,C0,epsA,T,f0,fQ)
% CONSTANTS
kB = 8.617e-5;              % Boltzmann constant kB = 8,617343? 10-5 eV/K
F = 96485.3;                % Faraday constant (C/mol)
eps0 = 8.854e-12;           % Vaccum permitivity [F/m]

% DERIVED PARAMETERS
% Inverse Debye's screening length [1/m]
kappa = sqrt(2*C0*F/(eps0*epsA*kB*T));
% Debye screening length [1/m]
Ld = 1./kappa;

% Average ion concentrations
cp_mean=0*r0;
cn_mean=0*r0;
n=length(r0); 
    for m=1:n
        if r0(m) < 100*Ld 
            % For small pore radii, use potential in cylindrical pore
            integrandp = @ (r) r.*exp(-f0*besselj(0,1i*r/Ld)/...
                (besselj(0,1i*r0(m)/Ld)*kB*T));
            integrandn = @ (r) r.*exp(f0*besselj(0,1i*r/Ld)/...
                (besselj(0,1i*r0(m)/Ld)*kB*T));
        else
            % For large pore radii, use approximation of plane surface
            integrandp = @ (r) r.*exp(-f0*exp((r-r0(m))/Ld)/(kB*T));
            integrandn = @ (r) r.*exp(f0*exp((r-r0(m))/Ld)/(kB*T));
        end
        intp= quadgk(integrandp,0,r0(m)); % Numerical integration
        intn= quadgk(integrandn,0,r0(m)); % Numerical integration
        % Average cation concentration [mol/m^3]
        cp_mean(m)=2*C0/(r0(m)^2)*intp;
        % Average anion concentration [mol/m^3] 
        cn_mean(m)=2*C0/(r0(m)^2)*intn; 
    end
% Consideration of the Stern layer based on the partition coefficient
% (Leroy et al., 2008). For fQ=0, no contribution of the Stern layer.
cp_mean=1/(1-fQ)*(cp_mean-C0*fQ);