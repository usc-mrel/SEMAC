function FourPlots_pixel(im_se,time,shown_img,colorange)
% Input im_se: [Nx Ny Nslice Ntime]
%       time: temporal vector
%       shown_img

figure;
subplot(2,3,[1 4]); imagesc(shown_img,colorange); axis image;
colormap(gray);
colorbar;
for idx = 1:1000
    [x,y,button] = ginput(1);
    row = floor(y)
    col = floor(x)
    subplot(2,3,2);
    plot(time, reshape(abs(im_se(row,col,:,:)), [length(time) 1]),'ro-');
%     hold on;
%     plot(time, ones(1,length(time)) * abs(I(row,col,11,11)),'bo-')
%     hold off
    %ylim([0 450])
    xline(0)
    ylabel('Magnitude','FontWeight','BOLD')
    xlabel('Freq (Hz)','FontWeight','BOLD')
    title('Magnitude','FontWeight','BOLD')
    set(gca,'XTick',time, 'XTickLabel',time)
    
    subplot(2,3,5);
    plot(time, reshape(angle(im_se(row,col,:,:))/pi*180, [length(time) 1]),'ro-');
    ylim([-180,180])
    xline(0)
    ylabel('Phase(degree)','FontWeight','BOLD')
    xlabel('Freq (Hz)','FontWeight','BOLD')
    title('Phase','FontWeight','BOLD')
    set(gca,'XTick',time, 'XTickLabel',time)
    
    subplot(2,3,3);
    plot(time, reshape(real(im_se(row,col,:,:)), [length(time) 1]),'ro-');
    xline(0)
    ylabel('Real','FontWeight','BOLD')
    xlabel('Freq (Hz)','FontWeight','BOLD')
    title('Real','FontWeight','BOLD')
    set(gca,'XTick',time, 'XTickLabel',time)
    
    subplot(2,3,6);
    plot(time, reshape(imag(im_se(row,col,:,:)), [length(time) 1]),'ro-');
    xline(0)
    ylabel('Imaginary','FontWeight','BOLD')
    xlabel('Freq (Hz)','FontWeight','BOLD')
    title('Imaginary','FontWeight','BOLD')
    set(gca,'XTick',time, 'XTickLabel',time)
    
    if button == 3
        return;
    end
end
end

