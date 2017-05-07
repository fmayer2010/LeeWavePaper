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


%% toggle plots
plotstream = 1;
plotvel = 1;
plotpres = 1;
plotflux = 1;
plotstability = 1;

%% initializations 

J = .3; % Juice parameter = Fr_outer^-1 = Fr_inner = u'/sqrt(g'*H)

U = 1; % undisturbed uniform horizontal velocity

N = 1; % undisturbed uniform horizontal buoyancy

H = J*U/N; % height of bathymetry ( = J when U=N=1)
Kc = N/U; % modulus of internal gravirty wave (IGW) wavenumber, = vertical  
% IGW wavenumber in hydrostatic limit, rad/L

Lc = 2*pi/Kc; % IGW (and vertical) wavelength

D = 3*Lc; % depth of plotted region 

L = 100*Lc;%0.01*2*pi*U./N; % wavelength of bathymetry, needs to be
% much larger than Lc for hydrostatic propagating solution to be valid

K = 2*pi/L; % wavenumber of bathymetry 

n = 100;

Lx = 2*L; % width of plotted domain

x = linspace(0,Lx,n);

z = linspace(-H,D,n);

[xx,zz] = meshgrid(x,z);


h = H.*cos(K.*x); % bathymetry
dhdx = -H*K.*sin(K.*x); % derivative of bathymetry

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
    -cos(K.*xx + Kc.*zz).^3);

%% generating stream functions
% psi = -U.*z_0 + eta(x,z_0) where eta(x,z_0) = delta(x,z)
psi_delta_0 = U.*(-zz + delta_0);
psi_delta_1 = U.*(-zz + delta_0 + J.*delta_1);
psi_delta_2 = U.*(-zz + delta_0 + J.*delta_1 + J^2.*delta_2);

psi_eta_0 = U.*(-zz + eta_0);
psi_eta_1 = U.*(-zz + eta_0 + J.*eta_1);
psi_eta_2 = U.*(-zz + eta_0 + J.*eta_1 + J^2.*eta_2);

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
    + 3.*cos(K.*xx+Kc.*zz).^2.*sin(K.*xx+Kc.*zz));

w1 = - H*K.*( sin(K.*xx+Kc.*zz) );
w2 = w1 + H*K*J.*( cos(2*K.*xx+Kc.*zz)...
    + sin(K.*xx+Kc.*zz).^2 - cos(K.*xx+Kc.*zz).^2);
w3 = w2 + H*K*J^2.*( -1/2.*sin(K.*xx+Kc.*zz)...
    - 3/2.*sin(K.*xx+Kc.*zz).*cos(2*K.*xx+Kc.*zz)...
    - 3/2.*cos(K.*xx+Kc.*zz).*sin(2*K.*xx+Kc.*zz)...
    + cos(2*K.*xx+2*Kc.*zz).*sin(K.*xx+Kc.*zz)...
    + 1/2.*sin(2*K.*xx+2*Kc.*zz).*cos(K.*xx+Kc.*zz)...
    + 3.*cos(K.*xx+Kc.*zz).^2.*sin(K.*xx+Kc.*zz));

%% velocity by diff

x_star = x./L;
z_star = z./Lc;

dx = x(end)-x(end-1);
dz = z(end)-z(end-1);
dx_star = dx/L;
dz_star = dz/Lc;

x1 = x_star;
z1 = linspace(z_star(1)+dz_star/2,z_star(end)-dz_star/2,n-1);
[xx1,zz1] = meshgrid(x1,z1);

x2 = linspace(x_star(1)+dx_star/2,x_star(end)-dx_star/2,n-1);
z2 = z_star;
[xx2,zz2] = meshgrid(x2,z2);

udiff = U.*(-1/dz.*diff(psi_eta_2,1,1));
wdiff = U.*(1/dx.*diff(psi_eta_2,1,2));


%% Pressure solution
rho = 1000;
% P0 = -1/2.*rho.*((u0-U).^2 + N^2.*eta_0.^2);
P1 = -rho*U.*(u1-U);
% higher order solution still mising some cross terms
P2 = -1/2.*rho.*(2*U.*(u2-U) + (u1-U).^2 + N^2.*(eta_0 + J.*eta_1).^2); 
P3 = -1/2.*rho.*(2*U.*(u3-U) + (u1-U).*(u2-U) + N^2.*(eta_0 + J.*eta_1).*(eta_0 + J.*eta_1 + J^2.*eta_2));

%% Energy flux via wave drag and form drag

WaveDrag1 = -u1.*w1;
WaveDrag2 = -u2.*w2;
WaveDrag3 = -u3.*w3;
Dwave1 = sum(WaveDrag1,2).*dx/2;
Dwave2 = sum(WaveDrag2,2).*dx/2;
Dwave3 = sum(WaveDrag3,2).*dx/2;

E1 = P1.*w1;
E2 = P2.*w2;
E3 = P3.*w3;
Ewave1 = sum(E1,2).*dx/2;
Ewave2 = sum(E2,2).*dx/2;
Ewave3 = sum(E3,2).*dx/2;



% Form drags
% D0 = -rho.*(1/2*u0.^2+N^2.^h.^2).*dhdx;
D1 = P1(1,:).*dhdx;
D1 = sum(D1).*dx/(x(end)/L)

D2 =  P2(1,:).*dhdx;
D2 = sum(D2).*dx/(x(end)/L)

D3 =  P3(1,:).*dhdx;
D3 = sum(D3).*dx/(x(end)/L)

Dstar = pi*rho*U^3/N*J^2

Dhonhydro  = Dstar*sqrt(1-U^2*K^2/N^2)

Eform = P3(1,:).*dhdx.*u3(1,:);
Eform = sum(Eform).*dx/(x(end)/L)
Estar = rho*U^4/N*pi*J^2




%% plotting stream functions
zz_star = zz./Lc;
xx_star = xx./L;
x_star = x./L;

if plotstream==1
    nlines = 6;
    psilines0 = [min(min(psi_eta_0)):(max(max(psi_eta_0))-min(min(psi_eta_0)))/nlines:max(max(psi_eta_0))];
    psilines1 = [min(min(psi_eta_1)):(max(max(psi_eta_1))-min(min(psi_eta_1)))/nlines:max(max(psi_eta_1))];
    psilines2 = [min(min(psi_eta_2)):(max(max(psi_eta_2))-min(min(psi_eta_2)))/nlines:max(max(psi_eta_2))];


    figure(1)
    subplot(3,1,1)
    % contour(xx_star,zz_star,psi_eta_0, 'LineColor','k', 'ShowText','on','TextStep',4);
    contour(xx_star,zz_star,psi_eta_0, 'LineColor','k','LevelList',psilines0);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title({['J=',num2str(J),', $\Psi_{\eta_0}$']},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    subplot(3,1,2)
    contour(xx_star,zz_star,psi_eta_1,'LineColor','k','LevelList',psilines1);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title('$\Psi_{\eta_1}$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    subplot(3,1,3)
    contour(xx_star,zz_star,psi_eta_2,'LineColor','k','LevelList',psilines2);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title('$\Psi_{\eta_2}$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')

    print('eta_solutions', '-depsc');

    psilines0 = [min(min(psi_delta_0)):(max(max(psi_delta_0))-min(min(psi_delta_0)))/nlines:max(max(psi_delta_0))];
    psilines1 = [min(min(psi_delta_1)):(max(max(psi_delta_1))-min(min(psi_delta_1)))/nlines:max(max(psi_delta_1))];
    psilines2 = [min(min(psi_delta_2)):(max(max(psi_delta_2))-min(min(psi_delta_2)))/nlines:max(max(psi_delta_2))];


    figure(2)
    subplot(3,1,1)
    contour(xx_star,zz_star,psi_delta_0, 'LineColor','k','LevelList',psilines0);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title({['J=',num2str(J),', $\Psi_{\delta_0}$']},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    subplot(3,1,2)
    contour(xx_star,zz_star,psi_delta_1, 'LineColor','k','LevelList',psilines1);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title('$\Psi_{\delta_1}$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    subplot(3,1,3)
    contour(xx_star,zz_star,psi_delta_2, 'LineColor','k','LevelList',psilines2);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'k')
    hold off
    title('$\Psi_{\delta_2}$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    print('delta_solutions', '-depsc');
end

%% Plotting velocities
if plotvel==1
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
    
end

%% Plotting Pressure

if plotpres==1
    figure(6)
    subplot(3,1,1)
    contourf(xx_star,zz_star,P1,'edgecolor','none');
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
    hold off
    title({['J=',num2str(J),', P to $O(J^1)$']},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,1,2)
    contourf(xx_star,zz_star,P2,'edgecolor','none');
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
    hold off
    title('P to $O(J^2)$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(3,1,3)
    contourf(xx_star,zz_star,P3,'edgecolor','none');
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
    hold off
    title('P to $O(J^3)$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
    
end

%% plotting Fluxes

if plotflux==1
    figure(7)
    subplot(2,1,1)
    plot(Dwave1./Dstar,z_star,Dwave2./Dstar,z_star,Dwave3./Dstar,z_star)
    title('Wave drag ($-uw(x,z)$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('$D/D_{linear}$','Interpreter','latex')
    legend('J^1','J^2','J^3')
    subplot(2,1,2)
    plot(Ewave1./Estar,z_star,Ewave2./Estar,z_star,Ewave3./Estar,z_star)
    title('Energy flux ($pw(x,z)$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('$E/E_{linear}$','Interpreter','latex')
    legend('J^1','J^2','J^3')


    figure(8)
    contourf(xx_star,zz_star,E2./Estar*Lx,'edgecolor','none');
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lc,'w')
    hold off
    title('E to $O(J^2)$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar
end



%% plotting maximum slope in delta as a function of J
Al = 0:0.01:1.1;
delta_max_0 = Al;
delta_max_1 = Al+1/2.*Al.^2;
delta_max_2 = Al+1/2.*Al.^2 + 1/2.*Al.^3;

Junstbl_delta = Al(find(delta_max_2>1,1))

dudzsq_2 = Al.^2;
dudzsq_3 = Al.^2 + Al.^3;
dudzsq_4 = Al.^2 + Al.^3 + 35/4.*Al.^4;
dudzsq_5 = Al.^2 + Al.^3 + 35/4.*Al.^4 - 9/2.*Al.^5;
dudzsq_6 = Al.^2 + Al.^3 + 35/4.*Al.^4 - 9/2.*Al.^5 + 81/4.*Al.^6;

Junstbl_RI = Al(find(dudzsq_6>4,1))

if plotstability==1
    figure(9)
    plot(Al,delta_max_0,Al,delta_max_1,Al,delta_max_2,Al,ones(size(Al)),'--k')
    title('maximum slope in $\delta$ as a function of J','Interpreter','latex')
    ylabel('$\partial_z\delta_{max}$','Interpreter','latex')
    xlabel('J','Interpreter','latex')
    legend('\delta_0','\delta_1','\delta_2','Location','nw')

    print('max_slopes', '-depsc');
    
    figure(10)
    plot(Al,dudzsq_2,Al,dudzsq_3,Al,dudzsq_4,Al,dudzsq_5,Al,dudzsq_6,Al,4.*ones(size(Al)),'--k')
    title('maximum shear squared $(\partial_z u)^2$ as a function of J','Interpreter','latex')
    ylabel('$(\partial_z u)^2_{max}$','Interpreter','latex')
    xlabel('J','Interpreter','latex')
    legend('dudz^2 to J^2','to J^3', 'to J^4','to J^5', 'to J^6','Location','nw')

    print('max_shearsq', '-depsc');
end