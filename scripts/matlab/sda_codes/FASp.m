function [f,U] = FASp(dt,xgtt)
% Nyquist frequency (highest frequency)
Ny = (1/dt)/2;
% number of points in xgtt
L = length(xgtt);
% Next power of 2 from length of xgtt
NFFT = 2^nextpow2(L);
% frequency spacing
df = 1/(NFFT*dt);
% Fourier amplitudes
U = abs(fft(xgtt,NFFT))*dt;
% Single sided Fourier amplitude spectrum
U = U(2:Ny/df+1);
% frequency range
f = linspace(df,Ny,Ny/df)';
end