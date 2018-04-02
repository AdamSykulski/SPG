load STAT.mat
load NSTAT.mat
MAT2=mean(1./(12*NSTAT(:,10:-1:1,:,3)),3)'; % average damping timescale (nonstationary method)
MAT1=mean(1./(12*STAT(:,10:-1:1,:,3)),3)'; % average damping timescale (stationary method)
%% FIGURE 5 IN PAPER
scrsz = get(0,'ScreenSize'); % The Figure size
Arthur5=figure('Position',[1 1 scrsz(3)/1 scrsz(3)/2.6]); subplot(1,2,1); imagesc(MAT1);            %# Create a colored plot of the matrix values
colormap(bone);  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
textStrings = num2str(MAT1(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:8,1:10);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(MAT1(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'YTick',1:10,...
        'YTickLabel',{'90','80','70','60','50','40','30','20','10','0'},...
        'XTick',1:8,...
        'XTickLabel',{'1','2','3','4','5','6','7','8'},...
        'TickLength',[0 0]);
    ylabel('v (cm/s)'); xlabel('1/\lambda (days)'); colorbar; title('(a)'); caxis([0.5 9.5]);
    subplot(1,2,2); imagesc(MAT2);            %# Create a colored plot of the matrix values
colormap(bone);  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
textStrings = num2str(MAT2(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:8,1:10);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(MAT2(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'YTick',1:10,...
        'YTickLabel',{'90','80','70','60','50','40','30','20','10','0'},...
        'XTick',1:8,...
        'XTickLabel',{'1','2','3','4','5','6','7','8'},...
        'TickLength',[0 0]);
    ylabel('v (cm/s)'); xlabel('1/\lambda (days)'); colorbar; title('(b)'); caxis([0.5 9.5]);
exportfig(Arthur5, 'Arthur5.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 13); % For exporting fig into paper