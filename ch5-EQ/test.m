g = zeros(18,1);
Qs = [0.5577,2.4277,1.707, 2.00636, 1.41075, 1.707,1.607,1.35183,1.407,1.97973, 2.0977, 1.43421, 1.207, 1.66776, 1.88196, 1.4267,  0.53541,  0.25372];
for i = 1:18
    g(i) = 5;
end
EQ(g,Qs)
