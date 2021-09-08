function dfc = dir_vec(f, tar_angle, angle, M, dis, c)
% 计算均匀线阵导向矢量
% Input:
%   f 频率 Hz
%   angle 目标方向 弧度
%   M 麦克风数量
%   dis 间距 m
%   c 声速 m/s
% Output:
%   dfc 导向矢量
 
dfc = zeros(M,1);
dt = dis / c;
cosd = cos(angle+pi/2-tar_angle);
for p = 1:M
    dfc(p) = exp(-1i*2*pi*f*dt*(p-1)*cosd);
end

end
