% Laprise and Peltier iterative algorithm to solve Long's Model

%% parameters
% 
N=1e-3; % upstream buoyancy frequency
U=1e-1; % upstream horizontal velocity
J=.8; % Juice number (internal froude number squared, N*h/U) 
h_0=J*U/N; % mximum height of bathymetry
kg = N/U; % lee wave wavenumber
Lg = 2*pi/kg; % lee wave wavelength

epsilon = 0.1; % launching number = nonhydrostatic number = U*kh/N
kh = epsilon*N/U; % wavenumber of hill
Lx = 2*pi/kh; % wavelength of hill

% some settings: 
%for witch with Nx = 53 and a grid 10*Lx wide:
% with epsilon=0.01: Jmax=0.81036 and deta_dz_max=0.6516
% with epsilon=0.1, Jmax=0.82079 and deta_dz_max=0.6619
% with epsilon=0.9, Jmax=1.16339 and deta_dz_max=0.9789
% with epsilon=1, jmax=1.2385 and deta_dz_max=1.0469

%% toggle plotting
% put a 1 in each element of the plot vector to get that plot
% plot vector: plotting=[psi, velocities, rho and N, P and fluxes, stability, system sketch]
% plotting=[1,1,1,1,1,1];
plotting=[1,0,0,0,0,1];
% plotting = [1,0,0,0,1,0];
% plotting=[1,1,0,0,0,0];

%% grid
% 
% % sine settings
% Nx = 151; % MUST BE ODD. # horiz cells. Also # sine components in synthesis.
% Nz = 150; % # vert cells in plots. No effect on streamline solution, but does effect derivatives.
% x = linspace(-2*Lx,2*Lx,Nx); % width of horiz domain. Make wider than hill!
% z = linspace(-h_0,Lg,Nz); % heigt of plotted domain.

% % witch settings
Nx = 513; % MUST BE ODD. # horiz cells. Also # sine components in synthesis.
Nz = 500; % # vert cells in plots. No effect on streamline solution, but does effect derivatives.
x = linspace(-15*Lx,15*Lx,Nx); % width of horiz domain. Make wider than hill!
z = linspace(0,1.5*Lg,Nz); % height of plotted domain.

[xx,zz]=meshgrid(x,z);
dx = x(2)-x(1);
dz = z(2)-z(1);
nyquist = 2*pi/2/dx; %maximum resolvable frequency with current grid
k = linspace(0,nyquist,(Nx+1)/2); % positive half of wave space
k = [k,-flip(k(2:end))]; % neg half of wave space (in the order of fft output)

% dimensionless grid
xx_star = xx./Lx;
zz_star = zz./Lg;
x_star = x./Lx;
z_star = z./Lg;
dx_star = dx/Lx;
dz_star = dz/Lg;

%% error tolerance
% Set when to consider solution converged
% Laprise and Peltier use 1 part in 1e6. 
tol = 1e-6*h_0;

%% bathymetry
% Sine bottom
% h = h_0.*cos(kh.*x);

% Witch of Agnesi (centered at x=0... so make sure x vector surrounds it)
h = h_0./(1+x.^2./Lx^2);


%% initial guess (linear solution)

eta_0 = h; % linearized bottom boundary condition
eta_hat = fft(eta_0); % fourier transform it

% % plot the eta_hat momentarily. This plot is written so it will be drawn over
% % later on... change the figure call to an integer>4 to keep the plot
% figure(1)
% plot(k,abs(eta_hat))
% title('$\hat{\eta}$','Interpreter','Latex')
% xlabel('k')
% ylabel('z')
% pause()

% To evaluate at z=h(x), must loop through each k, determine is evanescent
% or propagating, and then evalueate this k at each h(x). This generates an
% eta(k,h(x))... thus a matrix of size (Nx,Nx)
eta_hat_h = zeros(length(h),length(k)); 
eta_h = zeros(size(x));
for i = 1:length(k)
    if abs(k(i))<kg
        m=sign(U*k(i))*(kg^2-k(i)^2)^(1/2);
    else
        m=1i*(k(i)^2-kg^2)^(1/2);
    end
    
    for j = 1:length(h)
        eta_hat_h(j,i) = eta_hat(i).*(exp(1i.*( m.*h(j))));
    end
end

% now loop through the h's and ifft each one to get eta(x,h(x)). This will
% also be a matrix of size Nx*Nx, but we only need the diagonal since all
% the other terms are eta(x1,h(x2)) where x1~=x2.
for i = 1:length(h)
    eta_temp = real( ifft(eta_hat_h(i,:)));
    eta_h(i) = eta_temp(i);
end

% initial error
Error = eta_h-h;
Enorm = Error./h; % normalized by the height of the bathymetry

% % plot the error momentarily. This plot is written so it will drawn over
% % later on... change the figure call to an integer>4 to keep the plot
% figure(1)
% plot(x,Error)
% title('\eta error')
% xlabel('x')
% ylabel('z')
% pause()

max_error = max(abs(Enorm));

%% Iterative solution
% Update the bottom boundary condition with the old one minus Error, then
% repeat the fft->eval at h(x)->ifft->compute Error until Error is less
% than the tolerance.
% Stops after 1000 iterations to avoid infinite looping on a nonconv sol.
% Stops if new error > old error to avoid blowing up.
iter = 0;
old_error = 10*max_error;
while max_error>tol && iter<1000 && max_error<old_error
    
% update eta at the bottom with the error
    eta_0 = eta_0-Error;
% fft it
    eta_hat = fft(eta_0);
% evaluate at h(x)
for i = 1:length(k)
    if abs(k(i))<kg
        m=sign(U*k(i))*(kg^2-k(i)^2)^(1/2);
    else
        m=1i*(k(i)^2-kg^2)^(1/2);
    end
    
    for j = 1:length(h)
        eta_hat_h(j,i) = eta_hat(i).*(exp(1i.*( m.*h(j))));
    end
end

% ifft to eta(x,h(x))
for i = 1:length(h)
    eta_temp = real( ifft(eta_hat_h(i,:)));
    eta_h(i) = eta_temp(i);
end

% update Error, old_error, and iter
    Error = eta_h-h;
    Enorm = Error./h;
    old_error = max_error;
    max_error = max(abs(Enorm));
    iter=iter+1;
    if max_error>old_error
        disp(['WARNING: solution blowing up. Caught at iter=',...
            num2str(iter),', eta_0 last updated for iter=',...
            num2str(iter-1),', where max_error=',num2str(old_error)])
    end
end

%% Compute converged solution on grid
% Now that our eta_0 has converged to producing a good eta_h,
% compute eta_z (i.e., ifft(eta_hat(k,z evalujated on the z-grid)) )

% iterate through the k's to calculate the m's at each z
eta_hat_z = zeros(Nz,Nx);
for i = 1:length(k)
    if abs(k(i))<kg
        m=sign(U*k(i))*(kg^2-k(i)^2)^(1/2);
    else
        m=1i*(k(i)^2-kg^2)^(1/2);
    end
    
    for j = 1:length(z)
        eta_hat_z(j,i) = eta_hat(i).*(exp(1i.*( m.*z(j))));
    end
end

% ifft back to real space
eta_N = zeros(size(xx));
for i = 1:length(z)
    eta_N(i,:) = real( ifft(eta_hat_z(i,:)));
end


%% Stream function
% compute psi and plot it
psi = U.*(zz-eta_N);

% set how many streamlines to plot (actually plots one more than this 
% integer, and the bathymetry is the lowest streamline)
nlines = 3;
psilines = [U*h(1):(U*(max(z)-h(1)))/nlines:U*max(z)];

if plotting(1)==1
    figure(1)
    contour(xx_star,zz_star,psi, 'LineColor','k','LevelList',psilines);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./(2*pi*U/N),'k')
    hold off
    title({['$\Psi(x,z)$, J=',num2str(J),', $\epsilon$=',num2str(epsilon)]},...
        'Interpreter','latex')
    ylabel('z/$\delta$','Interpreter','latex')
    xlabel('x$k/(2\pi)$','Interpreter','latex')
end
    
%% Velocity
% compute u and w and plot them

% grids for differentials (that is, evaluaded half way between old grid
% locations
% for d/dz:
x1 = x_star;
z1 = linspace(z_star(1)+dz_star/2,z_star(end)-dz_star/2,Nz-1);
[xx1,zz1] = meshgrid(x1,z1);
% for d/dx:
x2 = linspace(x_star(1)+dx_star/2,x_star(end)-dx_star/2,Nx-1);
z2 = z_star;
[xx2,zz2] = meshgrid(x2,z2);

% compute u and z from derivatives of eta (something is buggy with 
% derivatives of psi... they are off by an order of magnitude
utot= U.*(1-1/dz.*diff(eta_N,1,1));
% utot= U.*(1/dz.*diff(psi,1,1));
uprime = U.*(-1/dz.*diff(eta_N,1,1));
w = U.*(1/dx.*diff(eta_N,1,2));

if plotting(2)==1
    figure(2)
    subplot(3,1,1)
    contourf(xx1,zz1,utot./(U),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$u_{total}/(U)$ via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
        subplot(3,1,2)
    contourf(xx1,zz1,uprime./(J*U),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$u_{prime}/(JU)$ via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,1,3)
    contourf(xx2,zz2,w./(epsilon*J*U),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$w/(\epsilon J U)$  via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
end

%% Density and Buoyancy
% compute rho and N_local from psi
g = 9.81;
rho0=1000;
% rho = rho0.*(1-N^2/U/g.*psi);
rho = rho0.*(1 - N^2/g.*(zz - eta_N));
rho_prime = rho-rho0;
rho_max = max(max(abs(rho)));
rho_scale = rho0*N^2*h_0/g;
% Nsq_local = N^2/U.*(1/dz.*diff(psi,1,1));
Nsq_local = N^2*(1-1/dz.*diff(eta_N,1,1));
N_local = sqrt(Nsq_local);

if plotting(3)==1
    figure(3)
    subplot(3,1,1)
    contourf(xx_star,zz_star,rho./rho0,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$\rho/\rho_0$ via $\Psi$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,1,2)
    contourf(xx_star,zz_star,rho_prime./rho_scale,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$\rho_{prime}/(\rho_0N^2h_0/g)$ via $\Psi$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,1,3)
    contourf(xx1,zz1,N_local./N,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('N(x,z)/N via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
end

%% Wave drag, pressure, and energy flux
P=rho0*U^2.*(1/dz.*diff(eta_N,1,1)).*(1+1/4.*(1/dz.*diff(eta_N,1,1)));

% interpolate wdiff to z1
wdiffint = (w(1:end-1,:)+w(2:end,:))./2;
% interpolate diff and P to x2
udiffint = (uprime(:,1:end-1)+uprime(:,2:end))./2;
Pint = (P(:,1:end-1)+P(:,2:end))./2;

wavedrag = rho0.*wdiffint.*udiffint;
energyflux = Pint.*wdiffint;

% integrate energy flux horizontally to get vertical energy flux
evert = sum(energyflux.*dx,2);

[xxint,zzint]=meshgrid(x2,z1);

if plotting(4)==1
    figure(4)
    subplot(3,2,[1,2])
    contourf(xxint,zzint,Pint./(rho0*U^2),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$p/(\rho U^2)$ via diff($\eta$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,2,[3,4])
    contourf(xxint,zzint,wavedrag./(rho0*pi*J^2*U^3/N/Lx),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('Dimensionless Wave drag = $\rho_0 uw/(\frac{\rho \pi J^2 U^3}{NL_{hill}}$) via diff($\eta$)',...
        'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,2,[5,6])
    contourf(xxint,zzint,energyflux./(rho0*pi*J^2*U^4/N/Lx),'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('Dimensionless Energy flux = $\rho_0 pw/(\frac{\rho \pi J^2 U^4}{NL_{hill}}$) via diff($\eta$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
%     subplot(3,2,6)
%     plot(evert,z1)
    
end

%% Stability
% convective instability
deta_dz = 1/dz.*diff(eta_N,1,1);
deta_dz_max = max(max(abs(deta_dz)))

% shear instability
du_dz = 1/dz.*diff(utot,1,1);
% interpolate N_local to du_dz
Nsq_int = (Nsq_local(1:end-1,:)+Nsq_local(2:end,:))./2;
% interpolate xx1 and zz1 to this too
xx3 = (xx1(1:end-1,:)+xx1(2:end,:))./2;
zz3 = (zz1(1:end-1,:)+zz1(2:end,:))./2;

Ri = Nsq_int./du_dz.^2;
Ri_min = min(min(Ri))
Rimax = Ri_min+10;
Ri(find(Ri>Rimax))=Rimax;


if plotting(5)==1
    figure(5)
        subplot(2,1,1)
    contourf(xx1,zz1,deta_dz,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title(['$\partial \eta/\partial z$ via diff($\Psi$), $\partial \eta/\partial z_{max}$=',...
        num2str(deta_dz_max)],'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
%     subplot(4,1,2)
%     contourf(xx3,zz3,N_int,'edgecolor','none')
%     hold on
%         fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
%     hold off
%     title('N via diff($\Psi$)','Interpreter','latex')
%     ylabel('z/Lc','Interpreter','latex')
%     xlabel('x/L','Interpreter','latex')
%     colorbar
%     subplot(4,1,3)
%     contourf(xx3,zz3,du_dz,'edgecolor','none')
%     hold on
%         fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
%     hold off
%     title('$\partial u/\partial z$ via diff(diff($\Psi$))','Interpreter','latex')
%     ylabel('z/Lc','Interpreter','latex')
%     xlabel('x/L','Interpreter','latex')
%     colorbar
    subplot(2,1,2)
    contourf(xx3,zz3,Ri,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title({['Richardson number via diff(diff($\Psi$)), not showing Ri$>$',num2str(Rimax)],...
        ['Ri$_{min}$=',num2str(Ri_min)]},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
end
    


%% plot for sketch of system
if plotting(6)==1
    figure('Color',[1,1,1])
    contour(xx_star,zz_star,psi, 'LineColor','k','LevelList',psilines);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./(2*pi*U/N),'k')
    hold off
    set(gca,'XTick','')
    set(gca,'YTick','')
    print('system_sketch', '-depsc');
end

%% ugly plots
% figure()
% plot(evert,z1)
% title('energy flux integrated horizontally')
% xlabel('flux')
% ylabel('z/Lg')