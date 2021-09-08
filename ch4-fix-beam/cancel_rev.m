clear; close all; clc;

% postion
c = 340;               
fs = 16000;
M = 8; % number of mic
r = [2.00 1.5 2; 
     2.03 1.5 2; 
     2.06 1.5 2; 
     2.09 1.5 2; 
     2.12 1.5 2; 
     2.15 1.5 2; 
     2.18 1.5 2; 
     2.21 1.5 2];         
s = [2.105 2 2];           
tar_angle = pi/2;

L = [5 4 3];              
beta = 0.69;
n = 4096;
h = rir_generator(c, fs, r, s, L, beta, n);

[x,fs] = audioread('speech.wav');
y = [];
for i = 1:M
    tmp = filter(h(i,:),[1.],x);
    filename = sprintf('speech_rev_%d.wav',i);
    audiowrite(filename, tmp, fs);
    y = [y ; tmp'];
end


% pi / 2
y1 = zeros(1,length(y(1,:)));
for i = 1:M
    y1 = y1 + y(i,:);
end

filename = 'speech_rev_ds_90.wav';
audiowrite(filename, y1, fs);
