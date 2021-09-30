% This file is a function that takes a time-series signal and returns the
% frequency and FFT (complex valued) of the signal.
%
% [f,Sig]=FFT_perso(t,sig)
%
% Last Modified: 11 October 2016

function [f,Sig]=FFT_perso(t,sig)

N=length(t);            % Number of points in signal
dt=t(2)-t(1);           % Sample time increment, (s) [assumes uniform sampling]
fs=1/dt;                % Sampling frequency, (Hz)

f=fs*(0:(N/2))/N;       % Postive part of frequency range, (Hz)
y_raw=fft(sig,N);       % Call Matlab FFT function
y_pos=y_raw(1:(N/2+1));   % Take the FFT values in postive frequencies
Sig=[y_pos(1);2*y_pos(2:end)]/N;  % Scale FFT values appropriately and output