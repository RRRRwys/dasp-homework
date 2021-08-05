c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [1 1.2 1];              % Receiver position [x y z] (m)
s = [1 1 1];                % Source position [x y z] (m)
L = [5 7 4];                % Room dimensions [x y z] (m)
beta = 0.69;                % Reverberation time (s)
n = 16384;                  % Number of samples


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
for i = [length(h)-2:-1:1]
    edc(i) = edc(i+1) + h(i) * h(i);
end
edc = 10 * log10(edc+eps);
figure(2);
plot(t,edc);
xlabel('Time(s)');
ylabel('Energy decay curve (dB)');