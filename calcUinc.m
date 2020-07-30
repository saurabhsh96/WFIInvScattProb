%Function to calculate the incident field x and y are meshgrid variables
%and src is the location of the src, kb is wavenumber of the background
%medium.
function uinc = calcUinc(x, y, src, kb)
    nu = 0;
    k = 2;
    x_dash = x-src(1);
    y_dash = y-src(2);
    rho = sqrt(abs(x_dash).^2 + abs(y_dash).^2);
    uinc = (-1j/4).*besselh(nu,k,kb.*rho);
end