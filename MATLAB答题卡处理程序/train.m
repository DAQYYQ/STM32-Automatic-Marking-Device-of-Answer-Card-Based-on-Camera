clear
clc
%%  训练字符
%%  读取样本数据
DATADIR='.\字母库\';                                   % 待处理图像目录 
dirinfo=dir(DATADIR);                                  % 获取图像目录所有文件信息
Name={dirinfo.name};                                   % 获取文件名
Name(1:2)=[];                                          % 去除文件夹固有信息
[nouse num_of_char]=size(Name);                        % 获取类别数量
count = 1;
images = [];
labels = [];
for  cnt=1  :num_of_char                               % for 循环读取所有文件夹
      cnt
      pathname=horzcat(DATADIR, Name{cnt},'\');        % 把路径和名字融合一起
      sub_dirinfo=dir(pathname);                       % 获取图像目录所有文件信息
      sub_Name={sub_dirinfo.name};                     % 获取文件名
      sub_Name(1:2)=[];  
      [nouse num_of_image]=size(sub_Name); 
      for i = 1: num_of_image
      image = imread(horzcat(pathname,sub_Name{i}));
%       灰度化
      if size(image,3) >1 
          image = rgb2gray(image);
      end
      bw  = im2bw(image,graythresh(image));
      %   最小化外接框 
      [bw2,BoundingBox1,im0] = edu_imgcrop3(bw,bw,image);
      bw = imresize(bw2,[28 28], 'bilinear');
%       把图像拉成向量
      bw1 = double(reshape(bw,28*28,1));
      images = [images,bw1];
%       字符对应的标签
      labels(count) = cnt;
      count = count +1;
      end
end
 
% d*n
%% 特征值归一化
[input,settings] = mapminmax(images);

 
%% 构造输出矩阵
s = length(labels) ;
output = zeros(s,num_of_char) ;
for i = 1 : s
   output(i,labels(i)) = 1;
end
output = output';

%% 设置神经网络参数并训练
% Create a Pattern Recognition Network
hiddenLayerSize = 200;
net = patternnet(hiddenLayerSize, 'trainscg');
%设置训练参数
net.trainparam.show = 50;
net.trainparam.epochs = 100 ;
net.trainparam.goal = 0.01 ;
net.trainParam.lr = 0.01 ;
 
%开始训练
%这里的 input 矩阵的行表示特征的维度，列代表一个样本；
%       output'每一列表示一个样本的标签。
[net,tr] = train(net,input,output) ;  
% 显示网络
% view(net)

y = net(input);
% 显示训练cost降低的过程
figure, plotperform(tr)
% 显示分类混淆矩阵
figure, plotconfusion(output,y) 
    
%% 保存网络模型
save model.mat settings net Name
%% 测试    
  Y = net(input);
[value,pred] = max(Y);
                    
aa = find(pred ==labels);
acc = length(aa)/length(labels)     
 
             
             