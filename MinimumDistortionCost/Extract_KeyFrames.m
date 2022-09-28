function [key_frames]=Extract_KeyFrames(image_histogram,numFrames,KeyFrame_T)
k=1;i=1;
frame_num=0;
key_frames(1)=i;
while(k<numFrames)
    gray_diff=image_histogram(k+i,:)-image_histogram(k,:);
    gray_sum=image_histogram(k+i,:)+image_histogram(k,:);
    pixel_diff1=zeros(1,256);
    for j=1:1:256
        if gray_sum(j)==0
             pixel_diff1(j)=0;
         else
             pixel_diff1(j)=(gray_diff(j)^2)/gray_sum(j);
        end
    end
%     pixel_diff1=(gray_diff.^2)./gray_sum;
    diff1=sum(pixel_diff1,2);
    if diff1>KeyFrame_T
        disp(['key-frame number=',num2str(k+i)]) 
        k=k+i;
        frame_num=frame_num+1;
        key_frames(frame_num+1)=k;
        i=1;
    else
        i=i+1;
    end
    if (k+i)>450
        break
    end
end


end
