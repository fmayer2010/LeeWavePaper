function [ h0 ] = spectralbathy( h_rms, k0, N, stdphi)
%returns a realization of Goff style spectral bathymetry using Zhao et al. 
% 2015's 1D recipe and values
% h_rms=root mean squared height in meters 
% k0 = large hill (corner) wavenumebr in cycles/meter (Zhao use off of 2500m)
% N = number of gridpoints in real space, make it odd

close all;

% discretization values
M = (N-1)/2; %number of gridpoints in frequency space (ideally a power of 2)

%topographic values
beta = 2.5; %power law
h_rms = 2*pi*h_rms; 
minwavelength = 312; % min IGW generating wavelength


% The range of relevant wave numbers should perhaps be determined by a
% maximum wavenumber using the nyquist wavenumber: 1/(2*dx) and a minimum
% wavenumber of 0. However, instead of the nyquist wave number, NF2014 use 
% the inverse of half the minimum IGW generating wavelength. I follow them
% here
kmax = 1/(minwavelength); %2*pi/(minwavelength); 
dx = 1/(kmax);%2*pi/(kmax)


%discretizing wavespace... 
k = linspace(0,kmax, M);


% defining the power spectrum of height (equals the fourier transform of h(x))
% equation from Zhao et al. 2015. Note that is has units of length
% squared...
PSD = h_rms^2 .* ( 1 + (k./k0).^2 ).^( -beta/2 );

% Convert power PSD to amplitude A(k),  A = SQRT[2 * PSD];
Amp = [PSD].^(1/2);

% make a random phase for each Amp component
random_phase = pi.* (4+stdphi.*randn(size(Amp)));%2*pi .* (rand(size(Amp)));% 

% Construct a frequency domain signal Z = Amp .* exp(i*random_phase)
Amp_rand = Amp .* exp(1i .* random_phase);

% make a two-sided Amp_rand that is conjugate symmetric (of size N now... 
% i.e. now we have it in negative to positive wavenumber space)
kfull = [flip(-1.*k), 0, k]; % full wavenumber spectrum
Amp_rand_sym = [conj(flip(Amp_rand)), abs(Amp_rand(1)), Amp_rand];

% iFFT (and fftshift) Amp_rand_sym for h
h = ifft(ifftshift(Amp_rand_sym)).*sqrt(N);
h0 = fftshift(h);


% %finally, this bathymetry is still so jagged. For now, I address this by
% %smoothing it directly:
% 
% figure()
% hh = Smooth(real(h0),3);
% plot(x,hh,'-.',x,real(h0))
% title(['h(x), dx = ',num2str(dx)])
% xlabel('horizontal location (meters)')
% ylabel('relative bathymetric height (meters)')
% legend('smoothed h(x)','h(x)')
% axis equal





end

