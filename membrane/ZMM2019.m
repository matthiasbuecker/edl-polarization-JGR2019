function s = ZMM2019(T,c0,mu,epsA,w,L1,R1,L2,R2,f0,fQ,wide)
    % ZMM2017 computes the effective conductivity of a sequence of two
    % cylindrical pores based on the 1D Marshall-Madden impedance (Marshall
    % and Madden, 1959) as proposed in Buecker et al. (2019)
    % Matthias Buecker, May 2019
    %
    % Input parameters:
    % T:        Temperature (K)
    % c0:       Bulk ion concentration (mol/m^3)
    % mu:       Ion mobility in bulk electrolyte (m^2/(Vs))
    % epsA:     Relative permitivity of electrolyte (-)
    % w:        Angular frequency (1/s)
    % L1:       Length of pore 1 (m)
    % R1:       Radius of pore 1 (m)
    % L2:       Length of pore 2 (m)
    % R2:       Radius of pore 2 (m)
    % f0:       Zeta potential (V) 
    % wide:     0 - no surface charge on walls of wide pore
    %           1 - surface charge on walls of wide pore

    % CONSTANTS
    kB = 8.617e-5;      % Boltzmann Constant (eV/K)
    F = 96485.3;        % Faraday Konstante (C/mol)
    eps0 = 8.854e-12;   % Vaccum permitivity (F/m)

    % MODEL PARAMETERS
    % Bulk conductivity (S/m)
    ka = 2*mu*c0*F;
    % Inverse Debye's screening length (1/m)
    kappa = sqrt(2*c0*F/(eps0*epsA*kB*T));
    
    % Surface charge densities in the diffuse layer (C/m^2)
    Sigmadp = 2*c0*F/kappa*(exp(-f0/(2*kB*T))-1);
    Sigmadm = 2*c0*F/kappa*(exp(f0/(2*kB*T))-1); 

    % Effective conductivities (S/m)
    kdp1 = 2*mu*Sigmadp/R1*wide;
    kdm1 = 2*mu*Sigmadm/R1*wide;
    kdp2 = 2*mu*Sigmadp/R2;
    kdm2 = 2*mu*Sigmadm/R2;
    % Mean ion concentrations zone 1 (mol/m^3) 
    cp1 = c0*(2*kdp1/ka+1);
    cn1 = c0*(2*kdm1/ka+1);
    % Mean ion concentrations zone 2 (mol/m^3)
    cp2 = c0*(2*kdp2/ka+1);
    cn2 = c0*(2*kdm2/ka+1);
    % Consideration of the Stern layer based on the partition coefficient
    % (Leroy et al., 2008). For fQ=0, no contribution of the Stern layer.
    cp1=1/(1-fQ*wide)*(cp1-c0*fQ*wide);
    cp2=1/(1-fQ)*(cp2-c0*fQ);

    % Total cross-sectional areas of cylindrical pores (m^2/pi)
    A1 = R1.^2;
    A2 = R2.^2;

    % Normalize mean ion concentrations and mobilities to largest
    % cross-sectional area (passive zone)
    Anorm = max(A1,A2);
    % (Effective) cation mobility in zone 1 (m^2/V/s)
    mp1=mu*cp1/c0.*A1./Anorm;
    % (Effective) anion mobility in zone 1 (m^2/V/s)
    mn1=mu*cn1/c0.*A1./Anorm;
    % (Effective) cation mobility in zone 2 (m^2/V/s)
    mp2=mu*cp2/c0.*A2./Anorm;
    % (Effective) anion mobility in zone 2 (m^2/V/s)
    mn2=mu*cn2/c0.*A2./Anorm;

    % Magnitude and phase of the Marshall-Madden impedance (Marshall and 
    % Madden, 1959)
    [absZ,phiZ]=ZMM(T,c0,w,L1,mp1,mn1,L2,mp2,mn2);
    
    % Convert 1D impedance to normalized effective conductivity (S/m)
    Zm = absZ.*(cos(phiZ)+1i*sin(phiZ));
    s = (L1+L2)./(Zm.*ka);