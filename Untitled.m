% 清理命令和变量，以免内存不够
clc;clear;

% 读数据
% 数据格式：x,y,z
x1 = load('betao.txt');
y1 = load('betat.txt');
z1 = load('loss.txt');

% 抽稀，以免内存不够
count    = 1;     % 新变量计数器
interval = 20; % 抽稀间隔
for i = 1 : interval : 30000
    x(count) = x1(i);
    y(count) = y1(i);
    z(count) = z1(i);
    count = count + 1;
end

%确定网格坐标（x和y方向的步长均取0.1）
[X,Y]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y)); 
%在网格点位置插值求Z，注意：不同的插值方法得到的曲线光滑度不同
Z=griddata(x,y,z,X,Y,'v4');
%绘制曲面
figure(1)
surf(X,Y,Z);
shading interp;
colormap(jet);
% view(0, 90);
colorbar;
title('epsilon=4')
for i=1:1000
    cla;
    hold on
    surf(X,Y,Z);
    shading interp;
    colormap(jet);
    plot3(x(i), y(i), z(i), 'ro-', 'Linewidth', 2);
    hold on
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if i==1
         imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,'DelayTime',1e-4);
    else
         imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',1e-4);
    end
end

%plot3(x, y, z, 'ro-', 'Linewidth', 2);
print(gcf, '-djpeg', 'xyz.jpg'); % save picture
