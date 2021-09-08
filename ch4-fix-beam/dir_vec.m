function dfc = dir_vec(f, tar_angle, angle, M, dis, c)
% �������������ʸ��
% Input:
%   f Ƶ�� Hz
%   angle Ŀ�귽�� ����
%   M ��˷�����
%   dis ��� m
%   c ���� m/s
% Output:
%   dfc ����ʸ��
 
dfc = zeros(M,1);
dt = dis / c;
cosd = cos(angle+pi/2-tar_angle);
for p = 1:M
    dfc(p) = exp(-1i*2*pi*f*dt*(p-1)*cosd);
end

end
