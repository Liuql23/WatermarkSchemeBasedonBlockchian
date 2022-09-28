function NC = NC(image_w,shuiyin)
%% 说明w1、w2为水印图像值
w1=double(image_w);
w2=double(shuiyin);
[sizew1,sizew2]=size(w1);
ee=w1;
sume=0;
for i=1:1:sizew1
    for j=1:1:sizew2
        ee(i,j)=w1(i,j)*w2(i,j);
        sume=sume+ee(i,j);
    end
end

sumw1=0;
for i=1:1:sizew1
    for j=1:1:sizew2
        sumw1=sumw1+w1(i,j)*w1(i,j);
    end
end
% sumw2=0;
% for i=1:1:sizew1
%     for j=1:1:sizew2
%         sumw2=sumw2+w2(i,j)*w2(i,j);
%     end
% end
NC=sume/sumw1;
end