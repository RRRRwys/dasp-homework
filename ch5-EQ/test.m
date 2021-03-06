clc;clear;clc;
g = zeros(18,1);
%Qs = [0.5577,2.4277,1.707, 2.00636, 1.41075, 1.707,1.607,1.35183,1.407,1.97973, 5.44574, 2.21293, 1.207, 1.66776, 0.88196, 2.4267,  1.53541,  1.25372];
Qs = [1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707,1.707];
% F0 = [55, 77, 110, 156, 220, 311,440,622,880,1200,
%       1800,2500,3500,5000,7000,10000,14000,20000];

for i = 1:18
    g(i) = -5;
end

%g(10) = 13.12668;
%g(7) = -30.19027;

%Qs(14) = Qs(14)/2;
%Qs(14) = Qs(14)*2;

[x,fs] = audioread('wnoise.wav');

y = EQ(g,Qs,fs,x);
y = y / max(abs(y));
audiowrite('eq_output.wav',y,fs);
