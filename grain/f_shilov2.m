function f = f_shilov2(omega,c0,D,e_m,e_p,a,zeta,T)
% F_SHILOV2 computes the coefficient gamma of the of the solution of the 
% model by Shilov and Dukhin (1970) as presented in Shilov et al. (2001).
%
% Matthias Buecker, May 2019
%
% Input paramters
% omega:        Angular frequency (rad/s) (scalar or array)
% c0:           Bulk ion concentration (mol/m^3)
% D:            Diffusion coefficient (m^2/s)
% e_m:          Relative permitivity of the medium (-)
% e_p:          Relative permitivity of the particle (-)
% a:            Particle radius (m)
% zeta:         Zeta potential (V)
% E0:           Strength of external electric field (V/m)
% T:            Absolute temperature (K)

% Fundamental constants
e_0 = 8.854e-12;                    % Vaccum permitivity [F/m]
F = 96485.3;                        % Faraday's constant [C/mol]
kB = 8.617e-5;                      % Boltzmann's constant [eV/K]

% Absolute dielectric constants (F/m)
e_m = e_m*e_0;
% Ion mobilities
mu = D/(kB*T);
% Bulk electrolyte conductivity (S/m)
ka = 2*mu*c0*F;
% Inverse Debye's screening length [1/m] [eq. (2)]
kappa = sqrt(ka/(D*e_m));
% Dimensionless potential [-] [eq. (28)]
V0 = zeta/(kB*T);
% Ionic surface conductivity of cations (S) [from eq. (27)]
kdp = ka/2/kappa*(2*(exp(-V0/2)-1));
% Ionic surface conductivity of anions (S) [from eq. (27)]
kdm = ka/2/kappa*(2*(exp(V0/2)-1));
% Dukhin number [eq. (33)]
Du = (kdp+kdm)/(ka*a);
% Relaxation time tau
t = a^2/(2*D);
% Relaxation time tau_alpha [eq. (35)]
ta = t*(Du+1)*ka^2/((2*kdp/a+ka)*(2*kdm/a+ka));

% Frequency-dependent coefficient
f = (2*Du-1)/(2*Du+2)-...
    3*ta/(2*t)*(2*kdp/a-2*kdm/a)^2/(ka*(2*Du+2))^2*...
    (1-1i*omega*ta./(1+sqrt(1i*omega*2*t)+1i*omega*ta));
end