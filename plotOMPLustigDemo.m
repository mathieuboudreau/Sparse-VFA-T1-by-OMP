function [] = plotOMPLustigDemo(rawllsT1Map, llsT1Map, maskPhantom, NRMSE_T1 )
%PLOTOMPLUSTIGDEMO Plots the resulting OMP reconstucted T1 maps for
%R=1,2,4,6 in a very nice figure. Not a flexible function.
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau

maskPhantom(:,:,2)=maskPhantom(:,:,1);
maskPhantom(:,:,3)=maskPhantom(:,:,1);

figure()
imagesc( cat(  2,cat(   2,rawllsT1Map.*maskPhantom(:,:,1)/max(max(rawllsT1Map)), llsT1Map(:,:,1).*maskPhantom(:,:,1)/max(max(max(llsT1Map.*maskPhantom)))  ),cat(2,llsT1Map(:,:,2).*maskPhantom(:,:,1)/max(max(max(llsT1Map.*maskPhantom))),llsT1Map(:,:,3).*maskPhantom(:,:,1)/max(max(max(llsT1Map.*maskPhantom))))));
axis image
colormap(jet)

text(0+195, 16,'R=1','Color','w', 'FontSize',26)
text(256+195, 16, 'R=2','Color','w', 'FontSize',26)
text(512+195, 16, 'R=4','Color','w', 'FontSize',26)
text(768+195, 16, 'R=6','Color','w', 'FontSize',26)

text(0+195, 240, cat(2,num2str(NRMSE_T1(4)*100,2),'%'),'Color','w', 'FontSize',26)
text(256+195, 240, cat(2,num2str(NRMSE_T1(1)*100,2),'%'),'Color','w', 'FontSize',26)
text(512+195, 240, cat(2,num2str(NRMSE_T1(2)*100,2),'%'),'Color','w', 'FontSize',26)
text(768+195, 240, cat(2,num2str(NRMSE_T1(3)*100,2),'%'),'Color','w', 'FontSize',26)

set(gcf,'outerposition',[1           1        1280         778])
end

