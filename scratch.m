dT = .1; %time step: 1/100*(period at cutoff frequency)
t = 0:dT:100;
T = 4;
w_o = 2*pi/T;
X = [sin(w_o*t)' sin(w_o*t)' sin(w_o*t)'];
d = 1;
figure
f_cutoff = 200*sqrt(w_o); %characteristic frequency from dispersion relation 

[FFT, WaveNum, Freq] = FFT2_grid_v2(X, d, dT);

% surf(2*pi*WaveNum, 2*pi*Freq, fliplr(abs(FFT)));
hold on
pcolor(2*pi*(WaveNum),2*pi*(Freq), flip(abs(FFT)));


xlim([-1.1*pi/d 1.1*pi/d])

xlabel('\mu')
ylabel('\Omega')
title('Unmodulated')
shading flat
Ccolormap('Seahawks')