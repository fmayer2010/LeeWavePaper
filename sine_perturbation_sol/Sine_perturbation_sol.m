%% streamlines over sine topography
% this script plots 0th, 1st, and 2dn order perterbation explansion in J =
% NH/U solutions to Long's Hydrostatic model of problem of lee waves over 
% the sinusoidal bathymetry h(x) = A*sin(k*x) as found in Smith 1977.

% that is, these solutions satisfy:
% (d/dz)^2*delta(x,z) + N^2/U^2*delta(x,z) = 0, 
% and
% delta(x,h(x)) = h(x) = delta_0(x,0) + ... + (J^n)*delta_n(x,0) + O(J)^(n+1), 
% where n is the order of the solution

% these solutions do not account for finite depth, although the plots have
% the parameter D for the sake of finite paper!

%% initializations 

J = .3; % Juice parameter = Fr_outer^-1 = Fr_inner = u'/sqrt(g'*H)

U = 1; % undisturbed uniform horizontal velocity

N = 1; % undisturbed uniform horizontal buoyancy

H = J*U/N; % height of bathymetry ( = J when U=N=1)
Kc = N/U; % modulus of internal gravirty wave (IGW) wavenumber, = vertical  
% IGW wavenumber in hydrostatic limit, rad/L

Lc = 2*pi/Kc; % IGW (and vertical) wavelength

D = 1*Lc; % depth of plotted region 

L = 1;%0.01*2*pi*U./N; % wavelength of bathymetry

K = 2*pi/L; % wavenumber of bathymetry 

N = 100;

x = linspace(0,2*L,N);

z = linspace(-H,D,N);

[xx,zz] = meshgrid(x,z);


h = H.*cos(K.*x); % bathymetry

%% oth order solution

delta_0 = H.*cos(K.*xx + Kc.*zz);

eta_0 = delta_0;

%% 1st order solution

delta_1 = 1/2*H.*sin(2*K.*xx + Kc.*zz);

eta_1 = delta_1 - 1/2.*H.*sin(2*K.*xx + 2*Kc.*zz);

%% 2nd order solution

delta_2 = +1/2*H.*cos(K.*xx + Kc.*zz);

eta_2 = delta_2 + H.*(1/2.* cos(K.*xx + Kc.*zz).*cos(2.*K.*xx + Kc.*zz) ...
    - 1/2.*sin(K.*xx + Kc.*zz).*sin(2.*K.*xx + Kc.*zz) ...
    + 1/2.*sin(K.*xx + Kc.*zz).*sin(2.*K.*xx + 2.*Kc.*zz) ...
    -cos(K.*xx + Kc.*zz).^2);

%% generating stream functions
% psi = -U.*z_0 + eta(x,z_0) where eta(x,z_0) = delta(x,z)
psi_delta_0 = -U.*zz + delta_0;
psi_delta_1 = -U.*zz + delta_0 + J.*delta_1;
psi_delta_2 = -U.*zz + delta_0 + J.*delta_1 + J^2.*delta_2;

psi_eta_0 = -U.*zz + eta_0;
psi_eta_1 = -U.*zz + eta_0 + J.*eta_1;
psi_eta_2 = -U.*zz + eta_0 + J.*eta_1 + J^2.*eta_2;

%% plotting stream functions
zz_star = zz./Lc;
xx_star = xx./L;
x_star = x./L;

figure(1)
subplot(3,1,1)
% contour(xx_star,zz_star,psi_eta_0, 'LineColor','k', 'ShowText','on','TextStep',4);
contour(xx_star,zz_star,psi_eta_0, 'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title({['J=',num2str(J),', $\Psi_{\eta_0}$']},'Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
subplot(3,1,2)
contour(xx_star,zz_star,psi_eta_1,'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title('$\Psi_{\eta_1}$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
subplot(3,1,3)
contour(xx_star,zz_star,psi_eta_2,'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title('$\Psi_{\eta_2}$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
xlabel('x/L','Interpreter','latex')

print('eta_solutions', '-depsc');

figure(2)
subplot(3,1,1)
contour(xx_star,zz_star,psi_delta_0, 'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title({['J=',num2str(J),', $\Psi_{\delta_0}$']},'Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
subplot(3,1,2)
contour(xx_star,zz_star,psi_delta_1, 'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title('$\Psi_{\delta_1}$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
subplot(3,1,3)
contour(xx_star,zz_star,psi_delta_2, 'LineColor','k');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
hold off
title('$\Psi_{\delta_2}$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
xlabel('x/L','Interpreter','latex')
print('delta_solutions', '-depsc');

%% velocity solutions

u0 = U;
u1 = u0 + J.*( sin(K.*xx+Kc.*zz) );
u2 = u1 - J^2.*( 1/2.*cos(2*K.*xx+Kc.*zz)...
    + sin(K.*xx+Kc.*zz).^2 - cos(K.*xx+Kc.*zz).^2);
u3 = u2 - J^3.*( -1/2.*sin(K.*xx+Kc.*zz)...
    - sin(K.*xx+Kc.*zz).*cos(2*K.*xx+Kc.*zz)...
    - cos(K.*xx+Kc.*zz).*sin(2*K.*xx+Kc.*zz)...
    + cos(2*K.*xx+2*Kc.*zz).*sin(K.*xx+Kc.*zz)...
    + 1/2.*sin(2*K.*xx+2*Kc.*zz).*cos(K.*xx+Kc.*zz)...
    + 2.*cos(K.*xx+Kc.*zz).*sin(K.*xx+Kc.*zz));

w1 = - H*K.*( sin(K.*xx+Kc.*zz) );
w2 = w1 + H*K*J.*( cos(2*K.*xx+Kc.*zz)...
    + sin(K.*xx+Kc.*zz).^2 - cos(K.*xx+Kc.*zz).^2);
w3 = w2 + H*K*J^2.*( -1/2.*sin(K.*xx+Kc.*zz)...
    - 3/2.*sin(K.*xx+Kc.*zz).*cos(2*K.*xx+Kc.*zz)...
    - 3/2.*cos(K.*xx+Kc.*zz).*sin(2*K.*xx+Kc.*zz)...
    + cos(2*K.*xx+2*Kc.*zz).*sin(K.*xx+Kc.*zz)...
    + 1/2.*sin(2*K.*xx+2*Kc.*zz).*cos(K.*xx+Kc.*zz)...
    + 2.*cos(K.*xx+Kc.*zz).*sin(K.*xx+Kc.*zz));

%% Plotting velocities

figure(3)
subplot(3,1,1)
contourf(xx_star,zz_star,u1,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title({['J=',num2str(J),', u to $O(J)$']},'Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
colorbar
subplot(3,1,2)
contourf(xx_star,zz_star,u2,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('u to $O(J^2)$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
colorbar
subplot(3,1,3)
contourf(xx_star,zz_star,u3,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('u to $O(J^3)$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
xlabel('x/L','Interpreter','latex')
colorbar

print('u_solutions', '-depsc');

figure(4)
subplot(3,1,1)
contourf(xx_star,zz_star,w1,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title({['J=',num2str(J),', w to $O(J)$']},'Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
colorbar
subplot(3,1,2)
contourf(xx_star,zz_star,w2,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('w to $O(J^2)$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
colorbar
subplot(3,1,3)
contourf(xx_star,zz_star,w3,'edgecolor','none');
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('w to $O(J^3)$','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
xlabel('x/L','Interpreter','latex')
colorbar

print('w_solutions', '-depsc');

%% velocity by diff

x_star = x./L;
z_star = z./Lc;

dx = x(end)-x(end-1);
dz = z(end)-z(end-1);
dx_star = dx/L;
dz_star = dz/Lc;

x1 = x_star;
z1 = linspace(z_star(1)+dz/2,z_star(end)-dz/2,N-1);
[xx1,zz1] = meshgrid(x1,z1);

x2 = linspace(x_star(1)+dx/2,x_star(end)-dx/2,N-1);
z2 = z_star;
[xx2,zz2] = meshgrid(x2,z2);

udiff = U.*(-1/dz.*diff(psi_eta_2,1,1));
wdiff = U.*(1/dx.*diff(psi_eta_2,1,2));

figure(5)
subplot(2,1,1)
contourf(xx1,zz1,udiff,'edgecolor','none')
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('u via diff($\Psi$)','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
set(gca,'XTickLabel','')
colorbar
subplot(2,1,2)
contourf(xx2,zz2,wdiff,'edgecolor','none')
hold on
    fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
hold off
title('w via diff($\Psi$)','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')
xlabel('x/L','Interpreter','latex')
colorbar

print('vel_via_diff', '-depsc');

%% plotting maximum slope in delta as a function of J
Al = 0:0.01:1;
delta_max_0 = Al;
delta_max_1 = Al+1/2.*Al.^2;
delta_max_2 = Al+1/2.*Al.^2 + 1/2.*Al.^3;

figure(6)
plot(Al,delta_max_0,Al,delta_max_1,Al,delta_max_2,Al,ones(size(Al)),'--k')
title('maximum slope in $\delta$ as a function of J','Interpreter','latex')
ylabel('$\partial_x\delta_{max}$','Interpreter','latex')
xlabel('J','Interpreter','latex')
legend('\delta_0','\delta_1','\delta_2','Location','nw')

print('max_slopes', '-depsc');

