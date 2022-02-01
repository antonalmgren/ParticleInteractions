function gif_spectrum(MM,m,r)

for i = 1:size(MM,3)
    NN(:,:,i) = MM(:,:,i)./m;
end


h = figure('Position', [300, 300, 500, 700]);
    set(gcf,'color','w');    
for i = 1:300%:size(NN,3)
    subplot(2,1,1)
    surface(squeeze(NN(:,:,i)));
    shading flat
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = '[# m^{-3}]';
    title(['Particle numbers, step ',num2str(i)])
    caxis([0 max(NN,[],'all')])
    set(gca,'ColorScale','log')
    %text(0, 1, labels(i),'Units','normalized')
    
    subplot(2,1,2)
    surface(squeeze(MM(:,:,i)));
    shading flat
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = '[\mu g m^{-3}]';
    title(['Particle numbers, step ',num2str(i)])
    caxis([0 max(MM,[],'all')])
    set(gca,'ColorScale','log')
    %text(0, 1, labels(i),'Units','normalized')    
    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'spectrum.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'spectrum.gif','gif','WriteMode','append','DelayTime',0.01);
    end
    clf
end
end