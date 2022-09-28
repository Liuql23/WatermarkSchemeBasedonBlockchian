function imshow_Filter_figure( lpdf,hpdf )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
x=0:1:15;




figure
subplot(211)
% scatter(x,hpdf,'k. ','linewidth',0.1);%bottom' | 'top' | 'origin'
stem(x,lpdf,'fill','.k');
set(gca,'xtick',[],'xticklabel',[])
set(gca,'Ytick',[-0.5:0.5:0.5])
set(gca,'XAxisLocation','origin')
title('h=Daubechies 8 wavelet decomp.low-pass','position',[8,-0.7],'fontname','Times New Roman','Color','k','FontSize',12)
box off 


subplot(212)
stem(x,hpdf,'fill','.k');
set(gca,'xtick',[],'xticklabel',[])
set(gca,'Ytick',[-0.5:0.5:0.5])
set(gca,'XAxisLocation','origin')
title('g=Daubechies 8 wavelet decomp.high-pass','position',[8,-1.3],'fontname','Times New Roman','Color','k','FontSize',12)
box off 


end

