clc; clear all; close all;
warning off all;
%% 程序功能：答题卡1识别
I = imread('photo_4.bmp');

%% 图像大小调整及其对比度增强 
if size(I, 1) > 2000
    I = imresize(I, 0.2, 'bilinear');
end
% 对比度增强 
I1 = imadjust(I, [0 0.6], [0 1]); 
figure 
subplot(2, 1, 1); imshow(I, []); title('原图像', 'FontWeight', 'Bold');
subplot(2, 1, 2); imshow(I1, []); title('对比度增强 ', 'FontWeight', 'Bold');
 
%% 图像滤波
% 高斯滤波窗口大小！
hsize = [3 3];
% 高斯滤波标准差大小 ！
sigma = 0.5; 
h = fspecial('gaussian', hsize, sigma);
I2 = imfilter(I1, h, 'replicate');
 
figure 
subplot(2,2, 1); imshow(I1, []); title('滤波前图像', 'FontWeight', 'Bold');
subplot(2,2, 2); imshow(I2, []); title('高斯后图像', 'FontWeight', 'Bold');
 % 中值滤波！
image1(:,:,1) = medfilt2(I1(:,:,1));   
image1(:,:,2) = medfilt2(I1(:,:,2));  
image1(:,:,3) = medfilt2(I1(:,:,3));  
subplot(223)
imshow(image1, []);
title('中值滤波') 

% 均值滤波！
a = [1 1 1                                
    1 1 1
    1 1 1];
a = a./9;
image2   = imfilter(I1, a, 'replicate');  
subplot(224)
imshow(image2, []);
title('均值滤波')
 
%%  灰度化
I3 = rgb2gray(I2);
figure;
subplot(221);imshow(I2);title('原图像') 
subplot(222);imshow(I3);title('灰度化') 

%% 二值化
% bw1 = im2bw(I3, graythresh(I3));
bw1 = im2bw(I3,0.62);
bw2 = ~bw1;
subplot(223);imshow(bw2);title('二值化') 

%% 边缘
edgebw = edge(I3);
subplot(224);imshow(edgebw);title('边缘图像') 

%% hough变换找直线   
[H, T, R] = hough(edgebw);
P = houghpeaks(H, 4, 'threshold', ceil(0.1*max(H(:))));
lines = houghlines(edgebw, T, R, P, 'FillGap', 20, 'MinLength', 200);
% 找最长的直线！
max_len = 0;
for k = 1 : length(lines)
    xy = [lines(k).point1; lines(k).point2]; 
    len = norm(lines(k).point1-lines(k).point2); 
    Len(k) = len;
    if len > max_len
        max_len = len;
        xy_long = xy;
    end
    XY{k} = xy; % 存储信息！
end
% 显示结果！
figure 
subplot(2, 2, 1); imshow(edgebw); title('边缘图像', 'FontWeight', 'Bold');
subplot(2, 2, 2); imshow(H, [], 'XData', T, 'YData', R, 'InitialMagnification', 'fit');
xlabel('\theta'); ylabel('\rho');
axis on; axis normal; title('霍夫变换域', 'FontWeight', 'Bold')
subplot(2, 2, 3); imshow(I1); title('原图像', 'FontWeight', 'Bold');
subplot(2, 2, 4); imshow(I1); title('最长直线标记', 'FontWeight', 'Bold');
hold on;
plot(xy_long(:,1), xy_long(:,2), 'LineWidth', 2, 'Color', 'b');
 
%%  根据直线计算倾斜角度
x1 = xy_long(:, 1);
y1 = xy_long(:, 2);
K1 = -(y1(2)-y1(1))/(x1(2)-x1(1));
angle = atan(K1)*180/pi;

%%  根据倾斜角度校正图像
I4 = imrotate(I1,   angle, 'bilinear');
bw3 = imrotate(bw2, angle, 'bilinear'); 
figure 
subplot(2, 2, 1); imshow(I1, []); title('原图像', 'FontWeight', 'Bold');
subplot(2, 2, 3); imshow(bw2, []); title('原二值图像', 'FontWeight', 'Bold');
subplot(2, 2, 2); imshow(I4, []); title('校正图像', 'FontWeight', 'Bold');
subplot(2, 2, 4); imshow(bw3, []); title('校正二值图像', 'FontWeight', 'Bold');
 
%%  形态学滤波
% 去除小面积区域！
bw4 = bwareaopen(bw3, round(0.01*numel(bw3)/100));
% 去除大面积区域 ！
bw4 = removelarge(bw4, round(0.035*numel(bw3)/100));

figure 
subplot(1, 2, 1); imshow(bw3, []); title('待操作图像', 'FontWeight', 'Bold');
subplot(1, 2, 2); imshow(bw4, []); title('滤波图像', 'FontWeight', 'Bold');
 
%% 连通域标记
pic2 = bw4;
[l,mm]=bwlabel(pic2,8);
index = [];
status =[];
test_set = [];
bound_srt = [];
m = 1;

%% 依次标记分割的字符bing 识别
load model.mat
figure
imshow(I4, []); title('识别结果', 'FontWeight', 'Bold');
results = [];
for i=1:mm 
%     分区域提取字符！
      temp = 0*pic2+1;
      [xxx,yyy] = find(l==i);
      for j = 1:length(xxx)
          temp(xxx(j),yyy(j)) = 0;
      end
%    膨胀后的最小化外接框 ！
     [bw2,BoundingBox1,im0] = edu_imgcrop3(temp,temp,I4);
%    最终的最小化外接框   ！
     BoundingBox = BoundingBox1;
%    面积大于 50才是字 且 白色占比不能太大。！
     if BoundingBox(3)*BoundingBox(4) > 50  && sum(sum(~bw2))/(BoundingBox(3)*BoundingBox(4)) <0.4 
%          imwrite(im0,['./字母库/' strcat(num2str(clock),'.jpg')])！
         status =[status;BoundingBox];
%          白底黑字！
%         figure;imshow(im0)
         rectangle('position',BoundingBox,'edgecolor','b');
          bw = imresize(bw2,[28 28], 'bilinear');
  %       把图像拉成向量！
         test_set = double(reshape(bw,28*28,1));
         %         神经网络预测函数predict！
         testInput = mapminmax('apply',test_set,settings);
         Y = net(testInput);
         [value,pred] = max(Y);
         text(BoundingBox(1),BoundingBox(2)-10, Name{pred}, 'color', 'b','fontsize', 12);
         results = [results ;Name{pred}];
     end
end
 
 daan = num2str(results)
 
  %% 计算得分
% 标准答案！
 biaozhundaan = {'A' 'A' 'B' 'C' 'D'};   
score = 0;     
count = strcmp(daan,biaozhundaan)
counts = sum(count)
% 每题分值为4分
score = counts * 4
disp(['总分：' num2str(score)])