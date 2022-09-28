function PSNR = PSNR(f1, f2)
%% 说明f1、f2为彩色图像
%% 函数
f1=double(f1);
f2=double(f2);
[size_f_1,size_f_2]=size(f1);

diff_f=f1-f2;
sum_diff=(sum(sum(diff_f(:,:,1).^2))+sum(sum(diff_f(:,:,1).^2))+sum(sum(diff_f(:,:,1).^2)))/(3*(size_f_1*size_f_2));
PSNR=20*log10(255/sqrt(sum_diff));

end