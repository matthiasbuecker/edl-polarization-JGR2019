% Visualization of spectra for membrane polarization
% Before you run this script, run parameters.m first
%
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
    [1,1,1]*0];     % Black

% Load constants and model parameters
load('./parameters','m')

% Zeta potential in wide pore
mstr = 'I';         % write 'J' for continuous or 'I' for discontinous EDL
if strcmp(mstr,'J')
    wide = 1;       % 1 - surface charge on walls of wide pore
else 
    wide = 0;       % 0 - no surface charge on walls of wide pore
end

% Loop over partitioning factors
for pp = [1,3:7]
    %% Analytical solution
    % High-resolution angular frequency array
    nf = 400;
    omeganf=logspace(-2,4,nf);
    % Surface charge density in the Stern layer (C/m^2)
    sigmaS = m.sigma*m.p(pp);
    % Surface charge density in the diffuse layer (C/m^2)
    sigmad = m.sigma-sigmaS;
    % Zeta potential (V)
    zeta = sigma2zeta(m.c0,m.epsA,sigmad,'d',m.T);
    % Inverse Debye's screening length (1/m)
    kappa = sqrt(2* m.c0*m.F/(m.eps0*m.epsA*m.kB*m.T));
    % Surface charge densities in the diffuse layer (C/m^2)
    Sigmadp = 2*m.c0*m.F/kappa*(exp(-zeta/(2*m.kB*m.T))-1);
    % Partition coefficient (-)
    fQ = sigmaS/(Sigmadp+sigmaS);
    % Original Buecker-Hoerdt model (Buecker and H?rdt, 2013)
    sMM1 = ZMM2013(m.T,m.c0,m.mu,m.epsA,omeganf,m.L1,m.R1,m.L2,m.R2,zeta,...
        fQ,wide);
    % The new approximation of the Buecker-Hoerdt model
    sMM2 = ZMM2019(m.T,m.c0,m.mu,m.epsA,omeganf,m.L1,m.R1,m.L2,m.R2,zeta,...
        fQ,wide);
    
    %% Numerical solution
    % Coupled polarization
    dat = dlmread(['.\model' mstr '\model' mstr '_spec0' num2str(pp-1) '.txt'],...
        '',5,0);
    s = dat(:,2);
    if strcmp(mstr,'J')
        datS = dlmread(['.\model' mstr '\model' mstr '_specS0' num2str(pp-1)...
            '.txt'],'',8,1);
        s = s+datS.';
    end
    omega = dat(:,1);    
    
    % Spline interpolation (for smoother eye-guiding lines)
    sreal = spline(omega,real(s),omeganf);
    simag = spline(omega,imag(s),omeganf);

    %% Visualization
    figure(1)
    % Real part
    semilogx(omeganf,sreal,'Color',c(1,:),'LineWidth',lw/3)
        hold on
    semilogx(omega,real(s),'Color',c(1,:),'Marker','.',...
        'MarkerSize',ms(1),'linestyle','none')
    if pp == 1
        semilogx(omeganf,real(sMM1),'LineWidth',lw*1.5,'Color','k',...
            'linestyle',':');
        semilogx(omeganf,real(sMM2),'LineWidth',lw*1.2,'Color','k',...
            'linestyle','-');
    end
    
    figure(2)
    % Imaginary part
    semilogx(omeganf,simag,'Color',c(1,:),'LineWidth',lw/3)
    hold on
    ph1 = semilogx(omega,imag(s),'Color',c(1,:),'Marker','.',...
        'MarkerSize',ms(1),'linestyle','none');
    if pp == 1
        ph2 = semilogx(omeganf,imag(sMM1),'LineWidth',lw*1.5,'Color','k',...
            'linestyle',':');
        ph3 = semilogx(omeganf,imag(sMM2),'LineWidth',lw*1.2,'Color','k',...
            'linestyle','-');
    end
end

figure(1)
xlim([1e-2,1e4])
% ylim([0.09,0.14]) % Model I
set(gca,'FontSize',fs,'XTick',10.^(-10:1:10))
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.03])
xlabel('\omega (rad/s)','FontSize',fs)
ylabel('\sigma''(\omega)/\sigma_a (-)','FontSize',fs)
pos = get(gca,'Position');
print('-dpng','-r600',['spectra_' mstr '_re'])

figure(2)
% set(lh,'Box','off') % Model I
xlim([1e-2,1e4])
set(gca,'FontSize',fs,'XTick',10.^(-10:1:10))
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.03])
ax = gca;
ax.YAxis.Exponent = -3;
xlabel('\omega (rad/s)','FontSize',fs)
ylabel('\sigma''''(\omega)/\sigma_a (-)','FontSize',fs)
set(gca,'Position',pos);
print('-dpng','-r600',['spectra_' mstr '_im'])