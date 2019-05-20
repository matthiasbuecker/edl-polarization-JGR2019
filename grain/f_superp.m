function f = f_superp(omega,c0,mu,muS,e_m,e_p,a,zeta,SigmaS,T)
% F_SUPERP computes the frequency-dependent reflection coefficient for the
% superposed Stern-layer, diffuse-layer, Maxwell-Wagner polarization process.
% The response is a superposition of the Schwarz (1962) and the Dukhin-
% Shilov (1974) responses.
%
% Matthias Buecker, May 2019
%
% Input paramters
% omega:        Angular frequency (rad/s) (scalar or array)
% c0:           Bulk ion concentration (mol/m^3)
% mu:           Ion mobility in bulk electrolyte (m^2/V/s)
% muS:          Ion mobility in Stern layer (m^2/V/s)
% e_m:          Relative permitivity of the medium (-)
% e_p:          Relative permitivity of the particle (-)
% a:            Particle radius (m)
% zeta:         Zeta potential (V)
% SigmaS:       Surface charge density of Stern layer (C/m^2)
% T:            Absolute temperature (K)

% Fundamental constants
e_0 = 8.854e-12;	% Vaccum permitivity (F/m)
F = 96485.3;      	% Faraday's constant (C/mol)
kB = 8.617e-5;    	% Boltzmann's constant (eV/K)

% Absolute dielectric constants (F/m)
e_m = e_m*e_0;
e_p = e_p*e_0;
% Diffusion coefficients [m^2/s]
D = mu*kB*T;
DS = muS*kB*T;
% Bulk electrolyte conductivity (S/m)
sigma0 = 2*mu*c0*F;
% Inverse Debye's screening length (1/m)
kappa = sqrt(sigma0/(D*e_m));
% Dimensionless potential (-)
V0 = zeta/(kB*T);
% Ionic surface conductivity of cations (S)
kdp = sigma0/2/kappa*(2*(exp(-V0/2)-1));
% Ionic surface conductivity of anions (S)
kdm = sigma0/2/kappa*(2*(exp(V0/2)-1));

% 1. Diffuse-layer polarization according to Dukhin and Shilov (1974)
% Dukhin number (-)
Du = (kdp+kdm)/(sigma0*a);
% Coefficient S (-)
S = (Du+1)*sigma0^2/((2*kdp/a+sigma0)*(2*kdm/a+sigma0));
% Relaxation time tau (s)
ta = a^2/(2*D)*S;
% Frequency-dependent coefficient (-)
fd = (2*Du-1)/(2*Du+2)-...
    3*S/(2)*(2*kdp/a-2*kdm/a)^2/(sigma0*(2*Du+2))^2*...
    (1-1i*omega*ta./(1+sqrt(1i*omega*2*ta/S)+1i*omega*ta));
% Complex conductivity of the diffuse layer (S)
sigmaD = sigma0*(1+2*fd)./(1-fd);

% 2. Stern-layer polarization according to Schwarz (1962)
% Surface conductivity of bound charges (S)
lambda0 = abs(SigmaS)*muS; 
% Time constant (s)
tS = a.^2./(2*DS);
% Complex conductivity of the Stern layer (S)
sigmaS = 1i*2*lambda0/a*omega*tS./(1+1i*omega.*tS);

% 3. Piecing everything together
% Effective complex conductivity of the sphere
sigmaC = sigmaD + sigmaS + 1i*omega*e_p;
% sigmaC = sigmaD + sigmaS + sigmaG;
% Complex conductivity of the electrolyte
sigmaA = sigma0 + 1i*omega*e_m;

% Dimensionless frequency-dependent reflection coefficient
f = (sigmaC-sigmaA)./(2*sigmaA+sigmaC);
end