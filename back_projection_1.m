% my back-projection reconstruction test

% 产生断层头模
I = phantom(512); % 产生512*512的头模
figure(1);
imshow(I, []);
title('512*512原始图像');

% 参数设置
[N, M]=size(I);
z = 2 * ceil(norm(size(I) - floor((size(I) - 1) / 2) - 1)) + 3; % radon变换默认平移点数/角度
                                                                % 直接打开matlab里面内置函数radon的文件
                                                                % 开头注释里面写有：如果调用时
                                                                % R =
                                                                % RADON(I,THETA,N)的N（在每个投影方向的平移数）没有指定
                                                                % 那么N=2*ceil(norm(size(I)-floor((size(I)-1)/2)-1))+3
Nt = 360; % 角度采样点数
Nd = N; % 平移数，实际上就是指，在theta方向上进行线积分的个数
x = pi / 180; % 角度增量
d = N / Nd; % 平移步长，即在固定角theta上移动的步长
theta = 1 : Nt;
a = zeros([N, M]);

% 求投影
[R, xp] = radon(I, theta); % 雷登变换以theta角度求投影

                           % If theta is a scalar, R is a column vector
                           % containing the Radon transform for theta
                           % degrees. 只是一个角度的投影
                           
                           % If theta is a vector, R is a matrix in which
                           % each column is the Radon transform for one of
                           % the angles in theta. 表示所有角度的线投影
                           
e = floor((z - Nd) / 2) + 2;
testR = R;
R = R(e:(Nd + e - 1), :); % R就是雷登变换的投影
%R1 = reshape(R, N, Nt);
figure(2);
imshow(R, []);

% R-L滤波器设计，离散的滤波器空间脉冲响应
g = -(Nd / 2 - 1):(Nd / 2); % g为[-255:256]
for i = 1:N
    if g(i) == 0
        hl(i) = 1 / (4 * d^2);
    else if  mod(g(i), 2) == 0 % 如果是偶数
            hl(i) = 0;
        else % 如果是奇数
            hl(i) = (-1) / (pi^2 * d^2 * (g(i)^2));
        end
    end
end
k = Nd / 2:(3 * Nd / 2 - 1); 

% 图像的卷积反投影重建
for m = 1:Nt % 360个角度迭代
 pm = R(:, m); % pm为投影数据
 u = conv(hl, pm); % 卷积滤波
    %pm = u(k); % 注释掉就可以进行直接反投影重建，不进行滤波
    Cm = ((N - 1) / 2) * (1 - cos((m - 1) * x) - sin((m - 1) * x));
    for i = 1:N
        for j = 1:N
            Xrm = Cm + (j - 1) * cos((m - 1) * x) + (i - 1) * sin((m - 1) * x);% 计算的是通过点(j - 1, i - 1) = (0, 0)的投影线所在的位置x'=Xrm
                                                                               % 之所以(j - 1)为xi，(i - 1)为yj，是因为j在内循环，以行循环进行的：(0, 0)->(1, 0)->(2, 0)->...
                                                                               % ((m - 1) * x)实际上就是角度phi
                                                                               % 计算Xrm就是要求pm的行号，累计投影值
                                                                               % 加上Cm的目的是进行原点平移，避免像素点的下标（index）出现负值
                                                                               % 原点最初是在图像的中心点(0, 0)
                                                                               % 平移后原点位于(N, 1)，即图像的左下角点(-(N+1)/2, -(N+1)/2)
            % 下面计算方法与Xrm + Cm一致
            %Xrm = (j - (N + 1) / 2) * cos((m - 1) * x) + (i - (N + 1) / 2) * sin((m - 1) * x) + (N + 1) / 2;
            if Xrm < 1 % 如果计算的投影所在x'的位置小于1（超出范围），直接令其为最小的位置x'=1
                n = 1;
                t = abs(Xrm) - floor(abs(Xrm));
            else % 如果计算的投影所在x'的位置大于或等于1（在x'轴范围以内）
                n = floor(Xrm);
                t = Xrm - floor(Xrm);
            end
            if n > (Nd - 1) % 如果计算的投影所在x'的位置大于(Nd - 1)（超出投影使计算的位置数），直接令其为最大的位置x'=Nd - 1
                n = Nd - 1;
            end
            p = (1 - t) * pm(n) + t * pm(n + 1);% 线性内插
                                                % p代表经过像素点(N + 1 - i, j)，即点(511, 0)的投影值
            a(N + 1 - i, j) = a(N + 1 - i, j) + p;% j在内循环，i在外循环，则为一行一行这样填充
        end
    end
end
figure(3);
imshow(a, []);

% recI = zeros([N, M]);
% for i = 1:N
%     for j = 1:N
%         for m = 1:Nt
%             recI(i, j) = recI(i, j) + testR(i * cos((m - 1) * x) + j * sin((m - 1) * x), m) * 1;
%         end
%     end
% end
% figure(3);
% imshow(recI, []);

% 用matlab内置的反投影直接重建
testI = iradon(testR, theta, 'linear','none');
figure(4);
imshow(testI, []);
