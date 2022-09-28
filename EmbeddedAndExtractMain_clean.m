clear
close all
clc


%% ��Ƶ����

% ������Ƶ���󣬲���Ϊ�ֲ�����
videofile='Suzie_cover.avi';
obj=VideoReader(videofile);

% ��ȡ��Ƶ֡����Ϣ
numFrames=obj.NumberOfFrames; %��Ƶ֡����
numHeight=obj.Height;%��Ƶ֡�ĳ���
numWidth=obj.Width;%��Ƶ֡�Ŀ��

%�������飬������Ƶ֡ͼ��ͻҶ�ֱ��ͼͳ����
image_data=uint8(zeros(numHeight,numWidth,3,numFrames)); %����numHeight*numWidth*3*numFrames��С��4-D����

gray_image_data=double(zeros(numHeight,numWidth,numFrames)); %����numHeight*numWidth*numFrames��С��3-D����

image_histogram=double(zeros(numFrames,256));
for k=1:numFrames
    image_data(:,:,:,k)=read(obj,k); %����k֡ͼ��洢��4-D������
    gray_image_data(:,:,k)=rgb2gray(image_data(:,:,:,k));%����k֡ͼ��ĻҶ�ͼ��洢��3-D������
    sequence=reshape(gray_image_data(:,:,k),numHeight*numWidth,1);
    [counts,centers]=hist(sequence,0:1:255,1);
    image_histogram(k,:)=counts;
end

imagew_data=image_data; % ˮӡ��Ƶ��֡ͼ�����ݾ���
gray_imagew_data=gray_image_data;% ˮӡ��Ƶ��֡�Ҷ�ͼ�����ݾ���
gray_test_imagew_data=gray_image_data;% ����ˮӡ��Ƶ��֡�Ҷ�ͼ�����ݾ���
%% ��Ƶ��������
rate=obj.FrameRate;% ��Ƶ֡����
scene_T=0.3; %��һ��ֱ��ͼ������ֵ,ȡֵΪ[0,1]
[numSceneFrameChange,norm_diff]=Scene_change_partitioning(image_histogram,numFrames,rate,scene_T); %��������
countScene=size(numSceneFrameChange,2)+1;

%% ��Ƶ�ؼ�֡��ȡ
KeyFrame_T=2400; %ֱ��ͼ������ֵ��һ��ȡֵ��[2000,inf)
[numKeyFrames]=Extract_KeyFrames(image_histogram,numFrames,KeyFrame_T); %��Ƶ�ؼ�֡��ȡ
Scene_KeyFrame=cell(1,countScene);

Scene_KeyFrame{1}=numKeyFrames(find(numKeyFrames < numSceneFrameChange(1)));
for i=2:1:countScene-1
    Scene_KeyFrame{i}=numKeyFrames(find(numKeyFrames >= numSceneFrameChange(i-1) & numKeyFrames < numSceneFrameChange(i)));    
end
Scene_KeyFrame{countScene}=numKeyFrames(find(numKeyFrames >= numSceneFrameChange(countScene-1)));

%% ��ȡˮӡ��Ϣ
w_lena_128=double(logical(imread('binary_lena_128.bmp')));
w_scene_128=double(logical(imread('binary_Scene_128.bmp')));
w_qust_64=double(logical(imread('binary_qust_64.bmp')));
w_mascot_64=double(logical(imread('binary_mascot_64.bmp')));
w_qust_32=double(logical(imread('binary_qust_32.bmp')));
w_qust_128=double(logical(imread('binary_qust_128.bmp')));

wsize0=32;wsize1=64;

w_lena_num1=w_lena_128(1:64,1:64);w_lena_num2=w_lena_128(1:64,65:128);
w_lena_num3=w_lena_128(65:128,1:64);w_lena_num4=w_lena_128(65:128,65:128);
w_scene_num1=w_scene_128(1:64,1:64);w_scene_num2=w_scene_128(1:64,65:128);
w_scene_num3=w_scene_128(65:128,1:64);w_scene_num4=w_scene_128(65:128,65:128);
W=zeros(128,128,10);
for Wi=1:10
%     W(:,:,Wi)=w_lena_128;
    W(:,:,Wi)=w_qust_128;
end
% W(:,:,1)=w_lena_num1;W(:,:,2)=w_lena_num2;   W(:,:,3)=w_lena_num3;W(:,:,4)=w_lena_num4;
% W(:,:,5)=w_scene_num1;W(:,:,6)=w_scene_num2; W(:,:,7)=w_scene_num3;W(:,:,8)=w_scene_num4;
% W(:,:,9)=w_qust_64;                          W(:,:,10)=w_mascot_64;
[m,n]=size(W(:,:,1));
% ALLW=[w_lena_128,w_scene_128,[w_qust_64;w_mascot_64]];
figure
subplot(141); imshow(w_lena_128);  title('lena','position',[60,154],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('lena')
subplot(142); imshow(w_scene_128); title('scene','position',[60,154],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('scene')
subplot(143); imshow(w_qust_64);   title('qust','position',[30,76],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('qust')
subplot(144); imshow(w_mascot_64); title('mascot','position',[30,76],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('mascot')

%% ˮӡ��Ϣ����
T=10;
Wm=W;
for i=1:10
    Wm(:,:,i)=Anorld(W(:,:,i),T,0);%W(:,:,i)Ϊˮӡͼ��TΪ������0ΪArnold���ܣ�1ΪArnold����
end

%% ���ܵ�ˮӡ��ϢSVD�ֽ�
for i=1:10
    [Uw(:,:,i),Sw(:,:,i),Vw(:,:,i)]=svd(Wm(:,:,i));
end


%% 2��D�����˲�������
% �����˲���������Get 2D wavelet filters - Daubechies 8
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768]; % 1D high pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);% 1D low pass decomposition filter
% construction of 2D wavelet filters
F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;
% �˲���ͼ����������ͼ��
imshow_Filter_figure( lpdf,hpdf )


%% ˮӡ��ϢǶ��

cost=cell(1,countScene);
result=cell(1,countScene);
index= cell(1,countScene);
aim_index=cell(1,countScene);
Uo=cell(1,countScene);
So=cell(1,countScene);
Vo=cell(1,countScene);

infoload=1.5; %% ��ϢǶ���غ� infoload >= 1.2ʱ�����ؼ�֡��ȡ��ˮӡ��ϢNCCֵ����1 
tic
for countScene_i=1:1:countScene

    current_scene=Scene_KeyFrame{countScene_i};
    countKeyFrame=size(current_scene,2); %���㵱ǰ�����йؼ�֡��
    
    sceneKeyFrameDWT=cell(countKeyFrame,4); sceneKeyFrameDWT_w=cell(countKeyFrame,4);
    for countKeyFrame_i=1:1:countKeyFrame        
        % ��ǰ�ؼ�֡R����С���ֽ�
        currentFrame=image_data(:,:,1,current_scene(countKeyFrame_i));
        [ca1,ch1,cv1,cd1]=dwt2(double(currentFrame),'haar');
        
        sceneKeyFrameDWT{countKeyFrame_i,1}=ca1;sceneKeyFrameDWT{countKeyFrame_i,2}=ch1;
        sceneKeyFrameDWT{countKeyFrame_i,3}=cv1;sceneKeyFrameDWT{countKeyFrame_i,4}=cd1;  
        
        % ����LH������Ӧ��ʧ���
        cost{1,countScene_i}{countKeyFrame_i}=DistortionCost(ch1);
        % Ѱ��Ŀ��С��ϵ��
        [ result{1,countScene_i}{countKeyFrame_i}, index{1,countScene_i}{countKeyFrame_i} ] = sort( cost{1,countScene_i}{countKeyFrame_i}(:)); 
        aim_index{1,countScene_i}{countKeyFrame_i}=index{1,countScene_i}{countKeyFrame_i}(1:m*n);
        aim_ch1=ch1(aim_index{1,countScene_i}{countKeyFrame_i});
        aim_ch1_matrix=reshape(aim_ch1,m,n);
        % Ŀ��С��ϵ��SVD�ֽ�
        [Uo{countScene_i}{countKeyFrame_i},So{countScene_i}{countKeyFrame_i},Vo{countScene_i}{countKeyFrame_i}]=svd(aim_ch1_matrix);
        % ˮӡ��ϢǶ��
        S=So{countScene_i}{countKeyFrame_i}+infoload.*Sw(:,:,countScene_i);
        aim_ch1_matrix_w=Uo{countScene_i}{countKeyFrame_i}*S*Vo{countScene_i}{countKeyFrame_i}';
        aim_ch1_w=reshape(aim_ch1_matrix_w,m*n,1);
        ch1_w=ch1;
        ch1_w(aim_index{1,countScene_i}{countKeyFrame_i})=aim_ch1_w;
        %С���ع�
        sceneKeyFrameDWT_w{countKeyFrame_i,1}=ca1;sceneKeyFrameDWT_w{countKeyFrame_i,2}=ch1_w;
        sceneKeyFrameDWT_w{countKeyFrame_i,3}=cv1;sceneKeyFrameDWT_w{countKeyFrame_i,4}=cd1;

        R=idwt2(double(ca1),double(ch1_w),double(cv1),double(cd1),'haar');
        imagew_data(:,:,1,current_scene(countKeyFrame_i))=uint8(R);
    end
end

%% ����PSNRֵ
Key_Frames_PSNR=zeros(1,numel(numKeyFrames));
for k=1:1:numel(numKeyFrames)
%     gray_imagew_data(:,:,numKeyFrames(k))=rgb2gray(image_data(:,:,:,numKeyFrames(k)));%����k֡ͼ��ĻҶ�ͼ��洢��3-D������
%     psnr=PSNR(gray_image_data(:,:,numKeyFrames(k)),gray_imagew_data(:,:,numKeyFrames(k))); 
    psnr=PSNR(image_data(:,:,:,numKeyFrames(k)),imagew_data(:,:,:,numKeyFrames(k))); 
    Key_Frames_PSNR(1,k)=psnr;
end
Key_Frames_PSNR_Mean=mean(Key_Frames_PSNR);%% ���йؼ�֡PSNRֵ��ֵ��




%% ˮӡ��Ϣ��ȡ
Wt_keyframes=cell(1,countScene);Wt=W;
Wt_temp=cell(1,countScene);
for countScene_i=1:1:countScene

    current_scene=Scene_KeyFrame{countScene_i};
    countKeyFrame=size(current_scene,2); %���㵱ǰ�����йؼ�֡��
    
    
    for countKeyFrame_i=1:1:countKeyFrame        
        % ��ǰ�ؼ�֡R����С���ֽ�
        currentFrame=imagew_data(:,:,1,current_scene(countKeyFrame_i));
        [ca1w,ch1w,cv1w,cd1w]=dwt2(double(currentFrame),'haar');
               

        % Ѱ��Ŀ��С��ϵ��
 
        aim_ch1w=ch1w(aim_index{1,countScene_i}{countKeyFrame_i});
        aim_ch1w_matrix=reshape(aim_ch1w,m,n);
        % Ŀ��С��ϵ��SVD�ֽ�
        [Ut,St,Vt]=svd(aim_ch1w_matrix);
        % ˮӡ��Ϣ��ȡ
        Swt=(St-So{countScene_i}{countKeyFrame_i})./infoload;
        Wmt=Uw(:,:,countScene_i)*Swt*Vw(:,:,countScene_i)';
        %ˮӡ��Ϣ���� 
        Wt_keyframes{countScene_i}{countKeyFrame_i}=Anorld(Wmt,T,1);
        binary_T=graythresh(Wt_keyframes{countScene_i}{countKeyFrame_i});
        Wt_keyframes{countScene_i}{countKeyFrame_i}=double(im2bw(Wt_keyframes{countScene_i}{countKeyFrame_i},binary_T));
    end
    
    for mi=1:m
        for nj=1:n
            for countKeyFrame_i=1:1:countKeyFrame
                Wt_temp{countScene_i}(countKeyFrame_i)=Wt_keyframes{countScene_i}{countKeyFrame_i}(mi,nj);
            end
            Wt(mi,nj,countScene_i)=mode(Wt_temp{countScene_i});
        end   
    end
  
end
toc
NC_value=zeros(1,10);
for NCi=1:10
    NC_value(NCi)=NC(W(:,:,NCi),Wt(:,:,NCi));
end
% w_lena_128_t(1:64,1:64)=Wt(:,:,1);w_lena_128_t(1:64,65:128)=Wt(:,:,2);
% w_lena_128_t(65:128,1:64)=Wt(:,:,3);w_lena_128_t(65:128,65:128)=Wt(:,:,4);
% 
% w_scene_128_t(1:64,1:64)=Wt(:,:,5);w_scene_128_t(1:64,65:128)=Wt(:,:,6);
% w_scene_128_t(65:128,1:64)=Wt(:,:,7);w_scene_128_t(65:128,65:128)=Wt(:,:,8);
% 
% w_qust_64_t=Wt(:,:,9);
% w_mascot_64_t=Wt(:,:,10);
% 
% figure
% subplot(141); imshow(w_lena_128_t);  title('lena','position',[60,154],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('lena')
% subplot(142); imshow(w_scene_128_t); title('scene','position',[60,154],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('scene')
% subplot(143); imshow(w_qust_64_t);   title('qust','position',[30,76],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('qust')
% subplot(144); imshow(w_mascot_64_t); title('mascot','position',[30,76],'fontname','Times New Roman','Color','k','FontSize',12)%gtext('mascot')
% NC_value(1)=NC(w_lena_128,w_lena_128_t);
% NC_value(2)=NC(w_scene_128,w_scene_128_t);
% NC_value(3)=NC(w_qust_64,w_qust_64_t);
% NC_value(4)=NC(w_mascot_64,w_mascot_64_t);

 figure
subplot(221)
imagesc(image_data(:,:,:,301));axis off;
subplot(222)
[LL,LH,HL,HH]=dwt2(image_data(:,:,1,301),'haar');
imagesc(LH);axis off;
subplot(223)
imagesc(cost{1,9}{1});colormap(gray(256));
coor=aim_index{1,9}{1}(:);x=rem(coor,180);y=(coor-x)/180+1;
hold on;plot(y,x,'b.','MarkerSize',0.1);hold off;axis off;
subplot(224)
imagesc(imagew_data(:,:,:,301));axis off;



