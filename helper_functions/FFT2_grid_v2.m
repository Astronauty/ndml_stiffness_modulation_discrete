function [FFT, W1, W2] = FFT2_grid_v2(data, dx, dy)

Nx = size(data, 2);
Ny = size(data, 1);


% nyFreqX = 1/2/dx;
% nyFreqY = 1/2/dy;
% 
% dFx = 1/(Nx*dx);
% dFy = 1/(Ny*dy);

if mod(Nx, 2) == 0
    freqX = 1/dx*(0:Nx-1)/Nx;
    freqX = freqX - freqX(Nx/2 + 1);
else
    freqX = 1/dx*(0:Nx-1)/Nx;
    freqX = freqX - freqX((Nx+1)/2);
end

if mod(Ny, 2) == 0
    freqY = 1/dy*(0:Ny-1)/Ny;
    freqY = freqY - freqY(Ny/2 + 1);
else
    freqY = 1/dy*(0:Ny-1)/Ny;
    freqY = freqY - freqY((Ny+1)/2);
end

% freqX = linspace(-nyFreqX, nyFreqX, Nx);
% freqY = linspace(-nyFreqY, nyFreqY, Ny);

% freqX = 1/dx*(0:Nx-1)/Nx;
% freqY = 1/dy*(0:Ny-1)/Ny;

% [W1, W2] = meshgrid(-nyFreqX:dFx:nyFreqX, -nyFreqY:dFy:nyFreqY);
[W1, W2] = meshgrid(freqX, freqY);
% W1 = fftshift(w1);
% W2 = fftshift(w2);
% FFT = 1/sqrt(Nx*Ny)*fliplr(fftshift(fft2(data)));
% FFT = 1/sqrt(Nx*Ny)*fliplr(fftshift(fft2(data)));

FFT = 1/sqrt(Nx*Ny)^2*fftshift(fft2(data));

