% ��������ͱ����������ڴ治��
clc;clear;

% ������
% ���ݸ�ʽ��x,y,z
x1 = load('betao.txt');
y1 = load('betat.txt');
z1 = load('loss.txt');

% ��ϡ�������ڴ治��
count    = 1;     % �±���������
interval = 20; % ��ϡ���
for i = 1 : interval : 30000
    x(count) = x1(i);
    y(count) = y1(i);
    z(count) = z1(i);
    count = count + 1;
end

%ȷ���������꣨x��y����Ĳ�����ȡ0.1��
[X,Y]=meshgrid(min(x):0.1:max(x),min(y):0.1:max(y)); 
%�������λ�ò�ֵ��Z��ע�⣺��ͬ�Ĳ�ֵ�����õ������߹⻬�Ȳ�ͬ
Z=griddata(x,y,z,X,Y,'v4');
%��������
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
