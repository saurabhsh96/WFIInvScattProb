%Function to create images of fields using vectors and fields
function createImage(x, y, uinc)
    figure();
    imagesc(x, y, real(uinc));  
    set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
    colorbar;
    axis equal tight;
    title('Real part of incident field on the object domain');
    xlabel('x','FontSize',13,'FontWeight','bold');
    ylabel('y','FontSize',13,'FontWeight','bold');

    figure();
    imagesc(x, y, imag(uinc));
    set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');    
    colorbar;
    axis equal tight;
    title('Imaginary part of incident field on the object domain');
    xlabel('x','FontSize',13,'FontWeight','bold');
    ylabel('y','FontSize',13,'FontWeight','bold');

    figure();
    imagesc(x, y, abs(uinc));
    set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
    colorbar;
    axis equal tight;
    title('Absolute value of incident field on the object domain');
    xlabel('x','FontSize',13,'FontWeight','bold');
    ylabel('y','FontSize',13,'FontWeight','bold');
end