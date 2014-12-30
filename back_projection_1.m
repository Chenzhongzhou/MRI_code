% my back-projection reconstruction test

% 产生断层头模
I = phantom(512); % 产生512*512的头模
figure(1);
imshow(I, []);
title('512*512原始图像');

% 参数设置
[N, M]=size(I);
z = 2 * ceil(norm(size(I) - floor((size(I) - 1) / 2) - 1)) + 3; % radon变换默认平移点数/角度
Nt = 360; % 角度采样点数
Nd = N; % 平移数
x = pi / 180; % 角度增量
d = N / Nd; % 平移步长
theta = 1 : Nt;
a = zeros([N, M]);

% 求投影
[R, xp] = radon(I, theta); % 雷登变换以theta角度求投影

                           % If theta is a scalar, R is a column vector
                           % containing the Radon transform for theta
                           % degrees. 只是一个角度一条线的投影
                           
                           % If theta is a vector, R is a matrix in which
                           % each column is the Radon transform for one of
                           % the angles in theta. 表示所有角度的线投影
                           
e = floor((z - Nd) / 2) + 2;
testR = R;
R = R(e:(Nd + e - 1), :); % R就是雷登变换的投影
R1 = reshape(R, N, 360);
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
    pm = u(k); % 注释掉就可以进行直接反投影重建，不进行滤波
    Cm = ((N - 1) / 2) * (1 - cos((m - 1) * x) - sin((m - 1) * x));
    for i = 1:N
        for j = 1:N
            Xrm = Cm + (j - 1) * cos((m - 1) * x) + (i - 1) * sin((m - 1) * x);
            if Xrm < 1
                n = 1;
                t = abs(Xrm) - floor(abs(Xrm));
            else
                n = floor(Xrm);
                t = Xrm - floor(Xrm);
            end
            if n > (Nd - 1)
                n = Nd - 1;
            end
            p = (1 - t) * pm(n) + t * pm(n + 1);
            a(N + 1 - i, j) = a(N + 1 - i, j) + p;
        end
    end
end
figure(3);
imshow(a, []);

% 用matlab内置的反投影直接重建
testI = iradon(testR, theta);
figure(4);
imshow(testI, []);
