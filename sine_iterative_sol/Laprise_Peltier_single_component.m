% Laprise and Peltier algorithm for a single sine wave bottom

%%
% parameters
N=1;
U=1;
J=0.55;
h_0=J*U/N;
kg = N/U;

epsilon = .01;
k = epsilon*N/U;

%%
% grid
Nx = 500;
Nz = 200;
x = linspace(0,4*pi/k,Nx);
z = linspace(-h_0,3*pi/kg,Nz);
[xx,zz]=meshgrid(x,z);
dx = x(2)-x(1);
dz = z(2)-z(1);

%%
% error tolerance
% Laprise and Peltier use 1 part in 1e6. 
tol = 1e-6;

%%
% bathymetry
h = h_0.*cos(k.*x);
% Lx = 2*pi/k;
% h = h_0./(1+x.^2./Lx^2);


%%
% initial guess
eta_0 = h;
eta_hat = fft(eta_0);

if k<kg
    m=sign(U*k)*(kg^2-k^2)^(1/2);
else
    m=1i*(k^2-kg^2)^(1/2);
end

eta_h = real( ifft(eta_hat .*(exp(1i.*( m.*h)))) );

% initial error

Error = eta_h-h;
Enorm = Error/h_0;


figure()
plot(x,Error)
title('\eta error')
xlabel('x')
ylabel('z')

max_error = max(abs(Enorm));
% compare tol, eval error, and Error

%%
% begin loop
iter = 0;
while max_error>tol && iter<1000
    
% update eta at the bottom with the error
    eta_0 = eta_0-Error;
    eta_hat = fft(eta_0);

% check if propagating or evanescent solution (this could be removed with
% the check above, but not too much overhead involved
if k<kg
    m=sign(U*k)*(kg^2-k^2)^(1/2);
else
    m=1i*(k^2-kg^2)^(1/2);
end

eta_h = real( ifft(eta_hat .*(exp(1i.*( m.*h)))) );
% eta_h = real( eta_0 .*(exp(1i.*(k.*x + m.*h))) );
% eta_h = eta_0 .*(cos(k.*x + m.*h));

% update error

    Error = eta_h-h;
    Enorm = Error/h_0;


max_error = max(abs(Enorm));

    iter=iter+1;
end

%%
% compute eta on the grid
eta_N = zeros(size(xx));
for i = 1:length(z)
    eta_N(i,:) = real( ifft(eta_hat .*(exp(1i.*( m.*z(i))))) );
end
% eta_N = real( (ones(length(z),1).*eta_0) .* (exp(1i.*(k.*xx + m.*zz))) );
% eta_N = (ones(length(z),1).*eta_0) .* (cos(k.*xx + m.*zz));

%%
% compute psi, u, w, rho, N, and plot them
psi = U.*(zz-eta_N);
zz_star = zz./(2*pi*U/N);
xx_star = xx./x(end);
x_star = x./x(end);

    nlines = 20;
    psilines = [min(min(psi)):(max(max(psi))-min(min(psi)))/nlines:max(max(psi))];

 
    figure(1)
    contour(xx_star,zz_star,psi, 'LineColor','k','LevelList',psilines);
    hold on
        fill([x_star(1),x_star, x_star(end)],[-h(1),h,-h(end)]./(2*pi*U/N),'k')
    hold off
    title({['$\Psi(x,z)$, J=',num2str(J)]},'Interpreter','latex')
    ylabel('z/Lc','Interpreter','latex')
    xlabel('x/L','Interpreter','latex')
    



