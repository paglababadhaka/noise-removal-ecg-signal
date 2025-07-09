clear all;
close all;
clc;

x = load('combine.mat');
w = x.val(1, 500:5000);
t = linspace(0,20,4501);

figure(1)
plot(w);
title('ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

%50 Hz eliminator
Fs = 500;                   % Sampling frequency
bw =5;                      % bandwidth

%Zero-calculation
Fr1 =50;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 50;                    % passband frequency
thetap = 2*pi*(Fp/Fs); r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num0 = poly([z1, z2]);
den0 = poly([p1, p2]);
[h0, f] = freqz(num0, den0, 4501, Fs);

figure(2);
subplot(211); 
plot(f, abs(h0)/max(abs(h0)), 'linewidth', 2);
grid on; 
title('Notch Filter 1');
xlabel('Frequency(Hz)');
ylabel('Magnitude Response');
subplot(212); 
plot(f, angle(h0)*180/pi, 'linewidth', 2);
grid on;
xlabel('Frequency(Hz)');
ylabel('Phase Response');

%100 hz eliminator
Fs = 500;            % Sampling frequency
Ts1 = 1/Fs;
bw =5;               % bandwidth

%Zero-calculation
Fr1 =100;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 100;                                     % passband frequency
thetap = 2*pi*(Fp/Fs); r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num1 = poly([z1, z2]);
den1 = poly([p1, p2]);
[h1, f] = freqz(num1, den1, 4501, Fs);

figure(3); 
subplot(211); 
plot(f, abs(h1)/max(abs(h1)), 'linewidth', 2);
grid on; 
title('Notch Filter 2');
xlabel('Frequency(Hz)');
ylabel('Magnitude Response');
subplot(212); 
plot(f, angle(h1)*180/pi, 'linewidth', 2); 
grid on;
xlabel('Frequency(Hz)'); 
ylabel('Phase Response');

%150 hz eliminator
Fs = 500;          % Sampling frequency
Ts1 = 1/Fs;
bw =5;             % bandwidth

%Zero-calculation
Fr1 =150;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); 
z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 150;              % passband frequency
thetap = 2*pi*(Fp/Fs); 
r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); 
p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num2 = poly([z1, z2]);
den2 = poly([p1, p2]);
[h2, f] = freqz(num2, den2, 4501, Fs);

figure(4); 
subplot(211);
plot(f, abs(h2)/max(abs(h2)), 'linewidth', 2);
grid on; 
title('Notch Filter 3');
xlabel('Frequency(Hz)');
ylabel('Magnitude Response');
subplot(212);
plot(f, angle(h2)*180/pi, 'linewidth', 2); 
grid on;
xlabel('Frequency(Hz)'); 
ylabel('Phase Response');

%Cascading
NUM1 = conv(num0, num1);
NUM = conv(NUM1, num2);

DEN1 = conv(den0, den1);
DEN = conv(DEN1, den2);

o = filter(NUM, DEN, w);

figure(5)
H = h0.*h1.*h2;
plot(f,H);
title('Cascaded Notch Filters');
xlabel('Frequency(Hz)'); 
ylabel('Magnitude Response');

figure(6)
plot(t, o);
title('Stage-1 Filtered Signal');
xlabel('Time'); 
ylabel('Magnitude');

%BAND PASS FILTER
fs = 500;

ecg_signal = o;

passband = [0.25 40];      % Passband frequencies in Hz
desired_attenuation = 60;  % Stopband attenuation in dB
beta = 5;                  % Kaiser window parameter, adjust as needed

%Design the filter
nyquist_freq = fs / 2;
normalized_passband = passband / nyquist_freq;

%Calculate the filter order using the Kaiser formula
delta_f = diff(normalized_passband);
A = -20*log10(10^(-desired_attenuation/20));
fir_order = ceil((A - 8) / (2.285*delta_f));

%Ensure the filter order is even for type 1 FIR filter
if rem(fir_order, 2) ~= 0
    fir_order = fir_order + 1;
end

%Generate the Kaiser window
kaiser_window = kaiser(fir_order+1, beta);

%Design the bandpass filter using the window method
n = -(fir_order/2):(fir_order/2);
h = (2 * passband(2) / fs) * sinc(2 * passband(2) * n / fs); % Ideal impulse response
h = h .* kaiser_window';                                     % Apply the Kaiser window

%Apply the filter to the ECG signal
filtered_ecg = conv(ecg_signal, h, 'same');

%Plot original and filtered signals
figure(7);
subplot(2,1,1);
plot(t, ecg_signal);
title('Original ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filtered_ecg);
title('Filtered ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

%Frequency Plotting
freq = -fs/2: fs/4501: (fs/2)-(fs/4501);
figure(8);
subplot(2,1,1);
plot(freq, abs(fftshift(fft(ecg_signal))));
title('Original ECG Signal');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0,50]);
subplot(2,1,2);
plot(freq, abs(fftshift(fft(filtered_ecg))));
title('Filtered ECG Signal');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0,50]);
sgtitle('Frequency Analysis')

%% For another value of row that is 30
x = load('combine.mat');
w = x.val(30, 500:5000);
t = linspace(0,20,4501);

figure(9)
plot(w);
title('ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

%50 Hz eliminator
Fs = 500;  % Sampling frequency
bw =5;     % bandwidth

%Zero-calculation
Fr1 =50;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 50;                                      % passband frequency
thetap = 2*pi*(Fp/Fs); r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num0 = poly([z1, z2]);
den0 = poly([p1, p2]);
[h0, f] = freqz(num0, den0, 4501, Fs);

%100 hz eliminator
Fs = 500;       % Sampling frequency
Ts1 = 1/Fs;
bw =5;          % bandwidth

%Zero-calculation
Fr1 =100;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 100;                                     % passband frequency
thetap = 2*pi*(Fp/Fs); r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num1 = poly([z1, z2]);
den1 = poly([p1, p2]);
[h1, f] = freqz(num1, den1, 4501, Fs);

%150 hz eliminator
Fs = 500;          % Sampling frequency
Ts1 = 1/Fs;
bw =5;             % bandwidth

%Zero-calculation
Fr1 =150;
theta1 = 2*pi*(Fr1/Fs);
z1 = exp(1j*theta1); z2 = exp(-1j*theta1);

%Pole-calculation
Fp = 150;                                     % passband frequency
thetap = 2*pi*(Fp/Fs); r = 1-(bw/Fs)*pi;
p1 = r*exp(1j*thetap); p2 = r*exp(-1j*thetap);

%Numerator vector & denominator vector
num2 = poly([z1, z2]);
den2 = poly([p1, p2]);
[h2, f] = freqz(num2, den2, 4501, Fs);

%Cascading
NUM1 = conv(num0, num1);
NUM = conv(NUM1, num2);

DEN1 = conv(den0, den1);
DEN = conv(DEN1, den2);

o = filter(NUM, DEN, w);

figure(10)
plot(t, o);
title('Stage-1 Filtered Signal');
xlabel('Time'); 
ylabel('Magnitude');

%BAND PASS FILTER
fs = 500;
ecg_signal = o;

passband = [0.25 40];      % Passband frequencies in Hz
stopband = [0 0.2 45 250]; % Stopband frequencies in Hz
desired_attenuation = 60;  % Stopband attenuation in dB
transition_width = 0.2;    % Transition width in Hz
beta = 5;                  % Kaiser window parameter, adjust as needed

%Design the filter
nyquist_freq = fs / 2;
normalized_passband = passband / nyquist_freq;
normalized_stopband = stopband / nyquist_freq;

%Calculate the filter order using the Kaiser formula
delta_f = diff(normalized_passband);
A = -20*log10(10^(-desired_attenuation/20));
fir_order = ceil((A - 8) / (2.285*delta_f));

%Ensure the filter order is even for type 1 FIR filter
if rem(fir_order, 2) ~= 0
    fir_order = fir_order + 1;
end

%Generate the Kaiser window
kaiser_window = kaiser(fir_order+1, beta);

%Design the bandpass filter using the window method
n = -(fir_order/2):(fir_order/2);
h = (2 * passband(2) / fs) * sinc(2 * passband(2) * n / fs); % Ideal impulse response
h = h .* kaiser_window'; % Apply the Kaiser window

%Apply the filter to the ECG signal
filtered_ecg = conv(ecg_signal, h, 'same');

%Plot original and filtered signals
figure(11);
subplot(2,1,1);
plot(t, ecg_signal);
title('Original ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, filtered_ecg);
title('Filtered ECG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

%Frequency Plotting
freq = -fs/2: fs/4501: (fs/2)-(fs/4501);
figure(12);
subplot(2,1,1);
plot(freq, abs(fftshift(fft(ecg_signal))));
title('Original ECG Signal');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0,50]);
subplot(2,1,2);
plot(freq, abs(fftshift(fft(filtered_ecg))));
title('Filtered ECG Signal');
xlabel('Frequency');
ylabel('Amplitude');
xlim([0,50]);
sgtitle('Frequency Analysis')