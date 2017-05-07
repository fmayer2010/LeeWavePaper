% computing the hydrostatic streamlines over a witch of agnesi (h(x)) using a fourier
% synthesis of component sine solutions, with each solved to second order 
% accuracy in internal froude number, J=Nh/U where H is the height of the
% kth component of h_hat(k) [H(k) = k(n)*h_hat(k(n))] and 
% h_hat is the fourier tranform of h(x)

J = 0.55; % background juice (note that each component of the solution is 
% accurate to its own juice, but this is the maximum of those, and thus
% represents the maximum error (error ~ J^(order+1)). It also determines
% the witch height.
U = 1; % undisturbed uniform horizontal velocity
N = 1; % undisturbed uniform horizontal buoyancy
hmax = J*U/N; % height of witch
Kc = N/U; % modulus of internal gravity wave (IGW) wavenumber, = vertical  
% IGW wavenumber in hydrostatic limit, rad/L
Lc = 2*pi/Kc; % IGW wavelength
Lw = 1*Lc; % half-length of witch. For hydrostatic flow, should be at least 
% 10 times larger than Lc
D = 3*Lc; % depth of plotted region... has no effect on solution
Nx = 300;
Nz = 40;
Lx = 10*Lw;
x = linspace(-Lx,Lx,Nx);
dx = x(2)-x(1);
z = linspace(0,D,Nz);
dz = z(2)-z(1);
[xx,zz] = meshgrid(x,z);
h = hmax./(1+((x)/Lw).^2); % witch bathymetry
dhdx = -hmax/dx.*(1./(1+((x+dx/2)/Lw).^2) - 1./(1+((x-dx/2)/Lw).^2)); % derivative of bathymetry

zz_star = zz./Lc;
xx_star = xx./Lw;
x_star = x./Lw;
z_star = z./Lc;



nyquist = 1/(2*dx);
Nk = 2*Nx+1;
% k = linspace(-nyquist,nyquist,Nk); % wave spectrum
% dk = k(2)-k(1);

% note: In Baines, the solution for the witch takes advantage of the
% evenness of the witch's transform and only integrates from 0 to infinity.
% Strangely, this seems to effect the result, as when I integrate from
% -infinity to infinity (or -nyquist to nyquist, really), I get a symmetric
% solutiuon over the witch, which is wrong. But when I change my k vector
% to go from 0 to nyquist, I get the right result.
k = linspace(0,nyquist,Nk+1); % wave spectrum
k = k(2:end); %get rid of the dc component
dk = k(2)-k(1);

h_hat = pi*hmax*Lw.*exp(-abs(k.*Lw)); % fourier transform of the witch 



%% test that fourier synthesis works
htest = zeros(size(x));
for i = 1:length(x)
    htest(i) = 1/(pi).*(sum(h_hat.*cos(k.*x(i)).*dk,2));
end

figure(3);
plot(x_star,h./Lc,x_star,htest./Lc,'--')
leg = legend('witch','synthesis of sine waves with $\hat{h}$');
set(leg,'Interpreter','latex')
title('testing fourier synthesis','interpreter','latex')
xlabel('x/Lw','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')

print('witch_bathymetry', '-depsc');


%% zeroth order solution
% testing fourier synthesis without the solution function, since the first
% order solution is so easy to write out

psitest = zeros(size(xx));
for i = 1:length(x)
    for j = 1:length(z)
        psitest(j,i) = -U*z(j)+U/(pi).*real(sum(h_hat.*...
            cos((k.*x(i) + Kc*z(j))).*dk));
    end

%     psitest(:,i) = U/(2*pi).*real(sum(( ...
%         (ones(Nz,1).*h_hat).*...
%         exp(1i.*( (ones(Nz,1).*k).*x(i) + Kc.*zz ))).*dk,2));
end


%% Solution using Nonhydrostatic Juice perturbation of component sine solutions

% Make to make a matrix of size [Nz,Nx,Nk] to fill with the component
% solution for each k
components = zeros(Nz,Nx,Nk);

for i = 1:length(k)
    for j = 1:length(x)
%         here I put in h_hat*dk as the component-wise height... I need a
%         component wise height for the Juice perturbations... using hmax
%         is wrong
    components(:,j,i) = (eta_nhs(k(i),N,U,h_hat(i)*dk,x(j),z,2))';
    end
end
    
psi_witch = zeros(size(xx));
for i = 1:length(x)
    for j = 1:length(z)
        psi_witch(j,i) = -U*z(j) + U/(pi).*sum(h_hat.*squeeze(components(j,i,:))'.*dk);
    end
end



%% Solution using hydrostatic Juice perturbation solution

% Make to make a matrix of size [Nz,Nx,Nk] to fill with the component
% solution for each k
components = zeros(Nz,Nx,Nk);

for i = 1:length(k)
    for j = 1:length(x)
    components(:,j,i) = (eta(k(i),N,U,h_hat(i)*dk,x(j),z,0))';
    end
end
    
psi_witch2 = zeros(size(xx));
for i = 1:length(x)
    for j = 1:length(z)
        psi_witch2(j,i) = -U*z(j) + U/(pi).*sum(h_hat.*squeeze(components(j,i,:))'.*dk);
    end
end



%% Plotting
nlines = D/Lc*4; %number of streamlines to plot

psilines = [min(min(psi_witch)):(max(max(psi_witch))-min(min(psi_witch)))/nlines:max(max(psi_witch))];

figure(4)
subplot(3,1,1)
contour(xx_star,zz_star,psi_witch, 'LineColor','k','LevelList',psilines);
hold on
    fill([x_star(1),x_star, x_star(end)],[min(z_star),h./Lc,min(z_star)],'k')
hold off
title({['$2^{nd}$ order non-hydrostatic fourier synthesis solution'],...
    ['J=',num2str(J),', $\epsilon$=',num2str(2*pi/Lw/Kc)]},'interpreter','latex')
set(gca,'XTickLabel','')
ylabel('z/Lc','Interpreter','latex')

psilines = [min(min(psi_witch2)):(max(max(psi_witch2))-min(min(psi_witch2)))/nlines:max(max(psi_witch2))];

subplot(3,1,2)
contour(xx_star,zz_star,psi_witch2, 'LineColor','k','LevelList',psilines);
hold on
    fill([x_star(1),x_star, x_star(end)],[min(z_star),h./Lc,min(z_star)],'k')
hold off
title({['$2^{nd}$ order hydrostatic']},'interpreter','latex')
set(gca,'XTickLabel','')
ylabel('z/Lc','Interpreter','latex')

psilines = [min(min(psitest)):(max(max(psitest))-min(min(psitest)))/nlines:max(max(psitest))];

subplot(3,1,3)
contour(xx_star,zz_star,psitest, 'LineColor','k','LevelList',psilines);
hold on
    fill([x_star(1),x_star, x_star(end)],[min(z_star),h./Lc,min(z_star)],'k')
hold off
title({'$0^{th}$ order hydrostatic and all propagating'},'interpreter','latex')
xlabel('x/Lw','Interpreter','latex')
ylabel('z/Lc','Interpreter','latex')

print('witch_eta', '-depsc');



%% Solution using delta form of Juice perturbation solution
% 
% % Make to make a matrix of size [Nz,Nx,Nk] to fill with the component
% % solution for each k
% components = zeros(Nz,Nx,Nk);
% 
% for i = 1:length(k)
%     for j = 1:length(x)
%     components(:,j,i) = (delta(k(i),N,U,h_hat(i)*dk,x(j),z,2))';
%     end
% end
%     
% psi_witch2 = zeros(size(xx));
% for i = 1:length(x)
%     for j = 1:length(z)
%         psi_witch2(j,i) = -U*z(j) + U/(pi).*sum(h_hat.*squeeze(components(j,i,:))'.*dk);
%     end
% end
% 
% psilines = [min(min(psi_witch2)):(max(max(psi_witch2))-min(min(psi_witch2)))/30:max(max(psi_witch2))];
% 
% figure()
% contour(xx_star,zz_star,psi_witch2, 'LineColor','k','LevelList',psilines);
% hold on
%     fill([x_star(1),x_star, x_star(end)],[min(z_star),h./Lc,min(z_star)],'k')
% hold off
% title({['$2^{nd}$ order hydrostatic fourier synthesis solution'],...
%     ['J=',num2str(J),', $\epsilon$=',num2str(2*pi/Lw/Kc),', using $\delta$']},interpreter','latex')
% xlabel('x/Lw','Interpreter','latex')
% ylabel('z/Lc','Interpreter','latex')
%         