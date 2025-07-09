clear all;
close all;
clc;

% Frequency components
f1 = 20 * (1 + 1);
f2 = 13 * (6 + 2);
f3 = 18 * (0 + 3);
f4 = 10 * abs(1 + 6 - 0 + 4);

% Time vector
fs = 1000; % Sampling frequency
t = 0:1/fs:2; % 2 seconds duration
N = length(t);
freq = -(fs/2): fs/N: (fs/2)-(fs/N);

% Generate signal with four frequency components
signal = cos(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*f3*t) + cos(2*pi*f4*t);

% Design a digital notch filter to remove the chosen frequency component
Fs = fs;
Ts1 = 1/Fs;
bw = 10;

% Zero-calculation
Fr = f2;
theta1 = 2*pi*(Fr/Fs);
z = exp(1j*theta1);
zp = exp(-1j*theta1);

% Pole-calculation
Fp1 = f2 - 5;
Fp2 = f2 + 5;
thetap1 = 2*pi*(Fp1/Fs);
thetap2 = 2*pi*(Fp2/Fs);
r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap1);
p1p = r*exp(-1j*thetap1);
p2 = r*exp(1j*thetap2);
p2p = r*exp(-1j*thetap2);
num = poly([z,zp]);
den = poly([p1, p1p, p2, p2p]);
[h, f] = freqz(num, den, 256, Fs);

% Apply the filter to the signal
filtered_signal = filter(num, den, signal);

% Plotting
% Time domain and Frequency domain plots of generated signal
figure(1);
subplot(2,1,1);
plot(t, signal);
title('Generated Signal (Time Domain)');

subplot(2,1,2);
plot(freq, abs(fftshift(fft(signal))));
title('Generated Signal (Frequency Domain)');

% Impulse response and frequency response of the designed filter
figure(2);
subplot(2,1,1);
impz(num, den);
title('Impulse Response of the Filter');

subplot(2,1,2);
freqz(num, den, fs);
title('Frequency Response of the Filter');

% Time domain and Frequency domain plots of the filtered signal
figure(3);
subplot(2,1,1);
plot(t, filtered_signal);
title('Filtered Signal (Time Domain)');

subplot(2,1,2);
plot(freq, abs(fftshift(fft(filtered_signal))));
title('Filtered Signal (Frequency Domain)');