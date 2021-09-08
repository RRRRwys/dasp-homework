clear; close all; clc;
c = 340;               
fs = 16000;              
r = [2 1.5 2];         
s = [2 2 2];           
L = [5 4 3];              
beta = 0.25;               
n = 4096;               
h = rir_generator(c, fs, r, s, L, beta, n);

%[sig,fs] = audioread('./res.wav'); 
%rs = fftfilt(h, sig);
%r2 = filter(h,1,sig);

t = [0:length(h)-1]./fs;
figure(1);
plot(t,h);
xlabel('Time(s)');
ylabel('Impulse response');

edc = zeros(1,length(h));
for i = [length(h)-1:-1:1]
    edc(i) = edc(i+1) + h(i) * h(i);
end
edc = 10 * log10(edc);
figure(2);
plot(t,edc);
xlabel('Time(s)');
ylabel('Energy decay curve (dB)');
