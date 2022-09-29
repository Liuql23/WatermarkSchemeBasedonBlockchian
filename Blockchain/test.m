clear all
close all
clc
% The size of watermark
W_H = 64;W_W = 64;

%╪сть cost
load cost                                                                                                               
cost_0 = cost{1,5}{1};
signe = zeros(size(cost_0));
% 
[C_H,C_W] = size(cost_0);
% [sort_cost,index] = sort(cost_0(:));
% aim_coord = index(1:W_H*W_W);
% signe(aim_coord)=1;
% imshow(cost_0,[])
% for iter1 = 1:C_H
%     for iter2 = 1:C_W
%         if signe(iter1,iter2)==1
%             hold on
%             plot(iter2,iter1,'b.','MarkerSize',1)
%         end
%     end
% end

load('.\Aim_coords\our_coord_5_1.mat')
signe(our_coord)=1;
imshow(cost_0,[])
for iter1 = 1:C_H
    for iter2 = 1:C_W
        if signe(iter1,iter2)==1
            hold on
            plot(iter2,iter1,'b.','MarkerSize',0.5)
        end
    end
end