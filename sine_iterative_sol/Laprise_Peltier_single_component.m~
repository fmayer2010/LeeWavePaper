% Laprise and Peltier iterative algorithm to solve Long's Model

%% parameters
% 
N=1e-3;
U=1e-1;
J=1;
h_0=J*U/N;
Lg = 
kg = N/U;
Lg = 2*pi/kg;

epsilon = 1.2;
kh = epsilon*N/U;

%% grid
% 
Nx = 53;
Nz = 60;
x = linspace(-10*pi/kh,10*pi/kh,Nx);
z = linspace(0,3*pi/kg,Nz);
[xx,zz]=meshgrid(x,z);
dx = x(2)-x(1);
dz = z(2)-z(1);
nyquist = 1/2/dx*2*pi;
k = linspace(0,nyquist,(Nx+1)/2);
k = [k,-flip(k(2:end))];

%% error tolerance
% 
% Laprise and Peltier use 1 part in 1e6. 
tol = 1e-6*h_0;

%% bathymetry
% 
Lx = 2*pi/kh;
% h = h_0.*cos(kh.*x);

h = h_0./(1+x.^2./Lx^2);


%% initial guess (linear solution)
% 
eta_0 = h;
eta_hat = fft(eta_0);
eta_hat_h = zeros(length(h),length(k));
eta_h = zeros(size(x));


% fill in the pisitive wavenumber results for m
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

% now fill in the second half of eta_hat_h with the reveresed first half
% minus the dc component
% eta_hat_h(:,length(k)+1:end) = eta_hat_h(:,2:length(k)-1);


% now loop through the h's and ifft each one to get eta(x,h(x))... this
% would be a matrix of size Nx*Nx. We only need the diagonal... since all
% the other terms are eta(x1,h(x2)) where x1~=x2.
for i = 1:length(h)
    eta_temp = real( ifft(eta_hat_h(i,:)));
    eta_h(i) = eta_temp(i);
end

% initial error

Error = eta_h-h;
Enorm = Error./eta_h;


% figure()
% plot(x,Error)
% title('\eta error')
% xlabel('x')
% ylabel('z')

max_error = max(abs(Error));
% compare tol, eval error, and Error

%% Iterative solution
%
iter = 0;
old_error = 10*max_error;
while max_error>tol && iter<1000 && max_error<old_error
    
% update eta at the bottom with the error
    eta_0 = eta_0-Error;
    eta_hat = fft(eta_0);

% check if propagating or evanescent solution (this could be removed with
% the check above, but not too much overhead involved
% if kh<kg
%     m=sign(U*kh)*(kg^2-kh^2)^(1/2);
% else
%     m=1i*(kh^2-kg^2)^(1/2);
% end
% 
% eta_h = real( ifft(eta_hat .*(exp(1i.*( m.*h)))) );
% % eta_h = real( eta_0 .*(exp(1i.*(k.*x + m.*h))) );
% % eta_h = eta_0 .*(cos(k.*x + m.*h));

% fill in the pisitive wavenumber results for m
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


% now loop through the h's and ifft each one to get eta(x,h(x))... this
% would be a matrix of size Nx*Nx. We only need the diagonal... since all
% the other terms are eta(x1,h(x2)) where x1~=x2.
for i = 1:length(h)
    eta_temp = real( ifft(eta_hat_h(i,:)));
    eta_h(i) = eta_temp(i);
end


% update error

    Error = eta_h-h;
    Enorm = Error./eta_h;

old_error = max_error;
max_error = max(abs(Error));

    iter=iter+1;
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


eta_N = zeros(size(xx));
for i = 1:length(z)
    eta_N(i,:) = real( ifft(eta_hat_z(i,:)));
end
% eta_N = real( (ones(length(z),1).*eta_0) .* (exp(1i.*(k.*xx + m.*zz))) );
% eta_N = (ones(length(z),1).*eta_0) .* (cos(k.*xx + m.*zz));

%% Stream function
% compute psi plot it
psi = U.*(zz-eta_N);
zz_star = zz./(2*pi*U/N);
xx_star = xx./Lx;
x_star = x./Lx;

    nlines = 19;
    psilines = [min(min(psi)):(max(max(psi))-min(min(psi)))/nlines:max(max(psi))];

 
    figure(1)
    contour(xx_star,zz_star,psi, 'LineColor','k','LevelList',psilines);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./(2*pi*U/N),'k')
    hold off
    title({['$\Psi(x,z)$, J=',num2str(J)]},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    
%% Velocity
% compute u and w and plot them
z_star = z./Lg;
dx_star = dx/Lx;
dz_star = dz/Lg;

x1 = x_star;
z1 = linspace(z_star(1)+dz_star/2,z_star(end)-dz_star/2,Nz-1);
[xx1,zz1] = meshgrid(x1,z1);

x2 = linspace(x_star(1)+dx_star/2,x_star(end)-dx_star/2,Nx-1);
z2 = z_star;
[xx2,zz2] = meshgrid(x2,z2);

udiff = U.*(1/dz.*diff(psi,1,1));
wdiff = U.*(-1/dx.*diff(psi,1,2));

figure(2)
    subplot(2,1,1)
    contourf(xx1,zz1,udiff./U,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('u/U via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(2,1,2)
    contourf(xx2,zz2,wdiff./U,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('w/U via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar


%% Density and Buoyancy
% compute rho and N_local from psi
g = 9.81;
rho = 1-N^2/U/g.*psi;
N_local = sqrt(N^2/U.*(1/dz.*diff(psi,1,1)));

figure(3)
    subplot(2,1,1)
    contourf(xx,zz,rho,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('$\rho/\rho_0$ via $\Psi$','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    set(gca,'XTickLabel','')
    colorbar
    subplot(2,1,2)
    contourf(xx1,zz1,N_local./N,'edgecolor','none')
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./Lg,'w')
    hold off
    title('N(x,z)/N via diff($\Psi$)','Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    colorbar

%% Momentum flux (wave drag)
% F = 