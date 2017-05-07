function [ Psi ] = componentpsi( K,N,U,h_hat,xx,zz,order )
%returns the dimensional vertical streamline of a leewave over the
%bathymetry h(x)=Hcos(kx)
%using the specified order accurate in J=Nh/U perturbation 
%expansion of the bottom boundary condition eta(xx,0)=h(x). This
%function is written such that it can be called in a fourier synthesis to
%produce streamfunctions for arbitrary topography

% the input variables are as follows:
% K is the wavelength of the sine function
% N is the constant upstream buoyancy frequency
% U is the constant upstream horizontal velocity
% h_hat is the Kth component of the fourier tranform of the bathymetry, in
% units of L^2
% xx and zz are matrices of the grid, produced via [xx,zz]=meshgrid(x,z)
% order is the order of accuracy of the solution. Can be 0, 1, or 2.

H = h_hat*K;
Kc = N/U;
J = N*H/U;

if K<=Kc %propagating solution

    eta0 = H.*cos(K.*xx + Kc.*zz);

    eta1 = 1/2*H.*sin(2*K.*xx + Kc.*zz) - 1/2.*H.*sin(2*K.*xx + 2*Kc.*zz);

    eta2 = 1/2*H.*cos(K.*xx + Kc.*zz)+...
        H.*(1/2.* cos(K.*xx + Kc.*zz).*cos(2.*K.*xx + Kc.*zz) ...
        - 1/2.*sin(K.*xx + Kc.*zz).*sin(2.*K.*xx + Kc.*zz) ...
        + 1/2.*sin(K.*xx + Kc.*zz).*sin(2.*K.*xx + 2.*Kc.*zz) ...
        -cos(K.*xx + Kc.*zz).^3);
else %evanescent solution (currently only solved to 0th order
    eta0 = H.*cos(K.*xx).*exp(-Kc.*zz);
    eta1 = 0;
    eta2 = 0;
end

if order==0
%     Psi = -U.*(zz + eta0);
    Psi = eta0;
elseif order==1
%     Psi = -U.*(zz +eta0 + J.*eta1);
    Psi = eta0 + J.*eta1;
elseif order==2
%     Psi = -U.*(zz +eta0 + J.*eta1 + J^2.*eta2);
    Psi = eta0 + J.*eta1 + J^2.*eta2;
else
    error('order must be 0, 1, or 2')
end





end

