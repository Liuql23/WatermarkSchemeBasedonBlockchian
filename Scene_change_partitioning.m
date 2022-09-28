function [scene_frames,norm_diff]=Scene_change_partitioning(image_histogram,numFrames,rate,threshold)

frames_diff=image_histogram(2:numFrames,:)-image_histogram(1:numFrames-1,:);
 current=image_histogram(1:numFrames-1,:);
 next=image_histogram(2:numFrames,:);
%  sum_bin=current+next;
 sum_bin=max(current,next);
 pixel_diff=frames_diff;
 for i=1:1:numFrames-1
     for j=1:1:256
         if sum_bin(i,j)==0
             pixel_diff(i,j)=0;
         else
             pixel_diff(i,j)=(frames_diff(i,j)^2)/sum_bin(i,j);
         end
     end
 end
diff=sum(pixel_diff,2);

%% 视频场景划分
norm_diff=(diff-min(diff))./(max(diff)-min(diff));
figure
bar(1:1:449,norm_diff,1)
xlabel('video frame number','FontName','Times New Roman','FontSize',12)
ylabel('$\hat{d}$', 'Interpreter', 'latex','FontName','Times New Roman','FontSize',16)
t=1;
scene=0;
T=ceil(rate*0.4);
for i=2:1:numFrames-2
    if norm_diff(i)>norm_diff(i-1) && norm_diff(i)>norm_diff(i+1)
        if norm_diff(i)>threshold
            t=t+1;
            jzd(t)=i;
            if jzd(t)-jzd(t-1)>T;
                scene=scene+1;
                disp(['scene changes frames=',num2str(jzd(t))]) 
                scene_frames(scene)=jzd(t);
            end
        end
    end
end
end
