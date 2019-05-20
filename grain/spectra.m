% Visualization of polarization spectra for spherical particles
% Before you run this script, run parameters.m first

% Matthias Buecker, May 2019
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear variables; close all; clc;

% Visualization settings
fs = 15;            % Font size
lw = 1;             % Line width
ms = [12,6,6];      % Marker size
c = [[1,1,1]*0.5;...% Grey 
    [1,1,1]*0.5;... % Grey 
    [1,1,1]*0];     % Black

% Load constants and model parameters
load('./parameters','m')

% Loop over partitioning factors
for pp = [1,3:length(m.p)]
    %% Analytical solution
    % Surface charge density in the Stern layer (C/m^2)
    sigmaS = m.sigma*m.p(pp);
    % Surface charge density in the diffuse layer (C/m^2)
    sigmad = m.sigma-sigmaS;
    % Zeta potential (V)
    zeta = sigma2zeta(m.c0,m.epsA,sigmad,'d',m.T);
    
    % High-resolution angular frequency array
    nf = 400;
    omega=logspace(-1,5,nf);
    % Frequency-dependent reflection coefficient (Shilov et al., 2001)
    fsh = f_shilov2(omega,m.c0,m.D,m.epsA,m.epsI,m.a,zeta,m.T);
    % Frequency-dependent reflection coefficient (Schurr, 1964)
    fsc = f_schurr(omega,m.c0,m.mu,m.muS,m.epsA,m.epsI,m.a,zeta,sigmaS,m.T);
    % Frequency-dependent reflection coefficient for combined
    % Schwarz-Dhukin-Shilov process
    fsp = f_superp(omega,m.c0,m.mu,m.muS,m.epsA,m.epsI,m.a,zeta,sigmaS,m.T);
    % Frequency-dependent reflection coefficient of Lyklema-Dhukin-Shilov
    % process
    fco = f_coupled(omega,m.c0,m.mu,m.muS,m.epsA,m.epsI,m.a,zeta,sigmaS,m.T);
    % Frequency-dependent reflection coefficient of coupled process
    % (empty Stern layer only)
    fcoD = f_coupled(omega,m.c0,m.mu,m.muS,m.epsA,m.epsI,m.a,zeta,sigmaS*0,m.T);
    % Effective normalized conductivity of the mixture (S/m)
  	scoD = (1+2*m.nu*fcoD)./(1-m.nu*fcoD);
    ssc = (1+2*m.nu*fsc)./(1-m.nu*fsc);
    ssp = (1+2*m.nu*fsp)./(1-m.nu*fsp);
    sco = (1+2*m.nu*fco)./(1-m.nu*fco);
    
    %% Numerical solution
    % Coupled polarization
    datG0 = dlmread(['.\modelG\modelG_spec00' num2str(pp-1) '.txt'],'',5,0);
    omegaG0 = datG0(:,1);
    sG0 = datG0(:,2);
    
    % Volumetric metal content of simulation
    L = 50e-6;
    nusim = 2*m.a.^3/(3*L.^3);
   
    % Convert to correct volumetric metal content
    fG0 = 1/nusim*(sG0-1)./(sG0+2);
    sG0 = (1+2*m.nu*fG0)./(1-m.nu*fG0);

    %% Visualization
    figure(1)
    % Real part
    semilogx(omega,real(ssp),'LineWidth',lw*1.5,'Color','k','linestyle',':');
    hold on
    plot(omega,real(sco),'LineWidth',lw,'Color','k','linestyle','-');
    plot(omegaG0,real(sG0),'Color',c(1,:),'Marker','.',...
        'MarkerSize',ms(1),'LineStyle','none')
    
    figure(2)
    % Imaginary part
    semilogx(omega,imag(ssp),'LineWidth',lw*1.5,'Color','k','linestyle',':');
    hold on
    plot(omega,imag(sco),'LineWidth',lw,'Color','k','linestyle','-');
    plot(omegaG0,imag(sG0),'Color',c(1,:),'Marker','.',...
        'MarkerSize',ms(1),'LineStyle','none')
    
    figure(3)
    % Imaginary part
    if pp == 1 || pp == 3 
    semilogx(omega,imag(scoD),'LineWidth',lw,'Color',c(pp,:),...
        'linestyle','--');
    hold on
    plot(omega,imag(sco),'LineWidth',lw,'Color',c(pp,:),'linestyle','-');
    plot(omega,imag(ssp),'LineWidth',lw,'Color',c(pp,:),'linestyle','-.');
    if pp == 3
        plot(omega,imag(ssc),'LineWidth',lw*1.5,'Color',c(pp,:),'linestyle',':');
    end
    plot(omegaG0,imag(sG0),'Color',c(pp,:),'Marker','.',...
        'MarkerSize',ms(1),'LineStyle','none','LineWidth',lw)
    end

end

figure(1)
xlim([1e-1,1e5])
set(gca,'FontSize',fs,'XTick',10.^(-10:1:10))
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.03])
xlabel('\omega (rad/s)','FontSize',fs)
ylabel('\sigma''(\omega)/\sigma_a (-)','FontSize',fs)
pos = get(gca,'Position');
print('-dpng','-r600','spectra_re')

figure(2)
lh = legend('Analytical (M=1)','Analytical (M>1)','Numerical',...
    'Location','Northeast');
set(lh,'Box','off')
xlim([1e-1,1e5])
set(gca,'FontSize',fs,'XTick',10.^(-10:1:10))
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.03])
xlabel('\omega (rad/s)','FontSize',fs)
ylabel('\sigma''''(\omega)/\sigma_a (-)','FontSize',fs)
set(gca,'Position',pos);
print('-dpng','-r600','spectra_im')

figure(3)
xlim([1e-1,3e4])
set(gca,'FontSize',fs,'XTick',10.^(-10:1:10))
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.03])
xlabel('\omega (rad/s)','FontSize',fs)
ylabel('\sigma''''(\omega)/\sigma_a (-)','FontSize',fs)
set(gca,'Position',pos);
print('-dpng','-r600','spectra_im_zoom')