clc;close all;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dir_vec

M = 8;
f = 8000;
dis = 0.03;
c = 340;
tar_angle = pi/2;

angles = 0: pi/1000 : pi*2.0;
g = zeros(length(angles),1);

for i = 1:length(angles)
    angle = angles(i);
    dirv = dir_vec(f, tar_angle, angle, M, dis, c);
    h = dirv;
    g(i) = abs(sum(h))/M;
end

figure
polar(angles', 20*log10(g) - min(20*log10(g)));
title('波束模式');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fq = 20:1000:20000;
wng = zeros(4,length(i));
for i = 1:length(fq)
    f = fq(i);
    d = dir_vec(f, tar_angle, tar_angle, 2, dis, c);
    h = d/2;
    wng(1,i) = abs(h'*d)^2/(h'*h);
    
    d = dir_vec(f, tar_angle, tar_angle, 4, dis, c);
    h = d/4;
    wng(2,i) = abs(h'*d)^2/(h'*h);
    
    d = dir_vec(f, tar_angle, tar_angle, 6, dis, c);
    h = d/6;
    wng(3,i) = abs(h'*d)^2/(h'*h);
    
    d = dir_vec(f, tar_angle, tar_angle, 8, dis, c);
    h = d/8;
    wng(4,i) = abs(h'*d)^2/(h'*h);
end

figure;
hold on;
plot(fq, 10*log10(wng(1,:)),'-r*');
hold on;
plot(fq, 10*log10(wng(2,:)),'-g*');
hold on;
plot(fq, 10*log10(wng(3,:)),'-bo');
hold on;
plot(fq, 10*log10(wng(4,:)),'-m+');
hold on;
ylabel('WNG(dB)');
xlabel('f(hz)')
title('白噪声增益');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 指向因子 DF

fq = 20:200:20000;
DF = zeros(4,length(i));
for i = 1:length(fq)
    f = fq(i);
    
    M = 2;
    gama = zeros(M);
    dt = dis/c;
    for ii = 1:M
        for jj = 1:M
            if(ii ~= jj)
                gama(ii,jj) = sin(2*pi*f*(jj-ii)*dt)/(2*pi*f*(jj-ii)*dt);
            else
                gama(ii,jj) = 1.0;
            end
        end
    end
    d = dir_vec(f, tar_angle, tar_angle, M, dis, c);
    h = d/M;
    DF(1,i) = abs(h'*d)^2/(h'*gama*h);
    
    M = 4;
    gama = zeros(M);
    dt = dis/c;
    for ii = 1:M
        for jj = 1:M
            if(ii ~= jj)
                gama(ii,jj) = sin(2*pi*f*(jj-ii)*dt)/(2*pi*f*(jj-ii)*dt);
            else
                gama(ii,jj) = 1.0;
            end
        end
    end
    d = dir_vec(f, tar_angle, tar_angle, M, dis, c);
    h = d/M;
    DF(2,i) = abs(h'*d)^2/(h'*gama*h);
    
    M = 6;
    gama = zeros(M);
    dt = dis/c;
    for ii = 1:M
        for jj = 1:M
            if(ii ~= jj)
                gama(ii,jj) = sin(2*pi*f*(jj-ii)*dt)/(2*pi*f*(jj-ii)*dt);
            else
                gama(ii,jj) = 1.0;
            end
        end
    end
    d = dir_vec(f, tar_angle, tar_angle, M, dis, c);
    h = d/M;
    DF(3,i) = abs(h'*d)^2/(h'*gama*h);
    
    M = 8;
    gama = zeros(M);
    dt = dis/c;
    for ii = 1:M
        for jj = 1:M
            if(ii ~= jj)
                gama(ii,jj) = sin(2*pi*f*(jj-ii)*dt)/(2*pi*f*(jj-ii)*dt);
            else
                gama(ii,jj) = 1.0;
            end
        end
    end
    d = dir_vec(f, tar_angle, tar_angle, M, dis, c);
    h = d/M;
    DF(4,i) = abs(h'*d)^2/(h'*gama*h);
end

figure;
hold on;
plot(fq, 10*log10(abs(DF(1,:))),'-r*');
hold on;
plot(fq, 10*log10(abs(DF(2,:))),'-g*');
hold on;
plot(fq, 10*log10(abs(DF(3,:))),'-bo');
hold on;
plot(fq, 10*log10(abs(DF(4,:))),'-m+');
hold on;
ylabel('DF(dB)');
xlabel('f(hz)')
title('指向因子');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
