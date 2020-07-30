%Final project WFI
clear;
close all;

%% Q1 and Q2 Learning about Hankel function and deriving usc
%Hankel is a third order Bessel Function given by hankle command in
%Matlab. Second question is done by substituting eq(1) and (2) in (3). 

%Taken from Matlab Documentation
k = 2;
nu = 0;
z = linspace(0.1,25,300);
H = besselh(nu,k,z);
figure(1);
plot(z,real(H),z,imag(H))
grid on
hold on
M = sqrt(real(H).^2 + imag(H).^2);
plot(z,M,'--')
legend('$J_0(z)$', '$Y_0(z)$', '$\sqrt{J_0^2 (z) + Y_0^2 (z)}$','interpreter','latex');
hold off;

%% Q3 Sketching the configuration

kb = 1;
lam = 2*pi/kb;
src = [lam/2, 10*lam];
domX = [0, 0, lam, lam, 0];
domY = [0, lam, lam, 0, 0];
figure(2);
line(domX, domY, 'LineWidth', 1.5); text(lam+0.1,lam,'(\lambda,\lambda)')
hold on;
scatter(src(1), src(2), 'filled'); text(lam/2+0.1,10*lam,'  \rho_s');
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
grid on;
title('Configuration of Source and Domain in Q3');
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');

%% Q5 Defining grid and the number of grid points
%Introducing grid
step = lam/20;
x_vec = 0:step:lam;
y_vec = 0:step:lam;
[x, y] = meshgrid(x_vec, y_vec);

%Numer of points
N = length(x).*length(y);

%Plotting grid
figure(2);
line(domX, domY, 'LineWidth', 1.5); text(lam+0.1,lam,'(\lambda,\lambda)')
line(x, y);
line(y, x);
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
title('Introducing a uniform grid in the object domain, h = \lambda/10');
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');

%% Q6 Calculation of incident field on the grid
uincInit = calcUinc(x, y, src, kb);
createImage(x_vec, y_vec, uincInit);

%% Q7 Variation in background field due to change in the location of the src and kb
%Changed position of source and changed kb
srcFar = [lam/2, 20*lam];
srcNear = [lam/2, 4*lam];
kbNew = 2;

%Field calculations - Far
% uincFar = calcUinc(x, y, srcFar, kb);
% createImage(x_vec, y_vec, uincFar);

%Field calculations - Near
% uincNear = calcUinc(x, y, srcNear, kb);
% createImage(x_vec, y_vec, uincNear);
 
% %Field calculations - Changed kb
uincKb = calcUinc(x, y, src, kbNew);
createImage(x_vec, y_vec, uincKb);

%% Q8 Introducting a contrast function of the object of my choice:
% I chose circle!! So, x^2 + y^2 = 0; Disc shape, however, due to the fact
% that the number of points in the grid are comparitively very few the
% circle won't be of a perfect shape.

center = [lam/2, lam/2];
expression = (x-center(1)).^2 + (y-center(2)).^2-(lam/6).^2;
%To make si >= 0
k_rho = ones(size(x)).*kb;
k_rho(expression <=0) = (1.2).*kb;
si_rho = (k_rho./kb).^2 - 1;

%% Q9 Making an image of the contrast function

figure(3);
imagesc(x_vec, y_vec, si_rho);
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
colorbar;
axis equal tight;
title('Contrast function of the object of choice');
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');

%% Q10 Defining the receiver domain L = 3*lam from (-lam, 1.5lam) to (2lam, 1.5lam)
%Adding the Rec domain to figure 1
drecX = [-lam 2*lam];
drecY = [1.5*lam 1.5*lam];

figure(2);
line(x, y);
line(y, x);
line(drecX, drecY, 'LineWidth', 1.5); hold on;
scatter(drecX, drecY, 'filled');
text(-lam-1.2,1.5*lam+1.8,'(-\lambda,1.5\lambda)');
text(2*lam-1.2,1.5*lam+1.8,'(2\lambda,1.5\lambda)');

title('Object Domain, Source and the Receiver Domain');

%% Q11 Introducing grid on Drec M > 0

Mdiv = 15; %Number of receivers
M = Mdiv + 1;
DrecX = -lam:3*lam/Mdiv:2*lam;
DrecY = 1.5*lam.*ones(size(DrecX));

figure(2);
scatter(DrecX, DrecY, 'LineWidth', 1.5);

%% Q12 Discretize data equation (appox. integral via sum)
%Mentioned in the report under question 12.

%% Q13 Reshape 2D array of contrast function; study reshape function

%Reshaped contrast function
siRs = reshape(si_rho, [N, 1]);

%% Q14 Building System Matrix A from discretized data equation
%Step size is same for x and y
delX = step;
delY = step;

%Covering the domain
nX = length(x_vec);
nY = length(y_vec);

%Defing the kth x
kX = 1:nX;
kY = 1:nY;
xK = (kX-1/2).*delX;
yK = (kY-1/2).*delY;

%Mid point meshgrid
[xK1, yK1] = meshgrid(xK, yK);

%Constant in front of sum in the discretized equation
const = -((kb^2)/16)*delX*delY;

%Finding system matrix 
A = zeros(M, N);
for ind = 1:M
    xm = DrecX(ind);
    ym = DrecY(ind);
    rhoMS = sqrt(abs(xm-xK1)^2 + abs(ym-yK1)^2);
    Greq = besselh(nu, k, (kb.*rhoMS));
    Greq = Greq*(calcUinc(xK1, yK1, src, kb).*(-4./1j));
    Greq = reshape(Greq, [1, N]);
    A(ind, :) = const.*Greq;
end

%% Q15 Compute SVD of system matrix 
%plot singular values, change M to invesigate dependance of singular values
%on M

%SVD
s = svd(A);
figure(4);
plot((s./max(s)), 'LineWidth', 1.5);
title(['Normalized Singular Values of the System Matrix, for M= ', num2str(M)]);
grid on;
%s > 0 -> eigen values are non-negative -> solutions available? Read more
%on this.. 

%% Q16 Computing scattered field

usc = A*siRs;

%% Q17 Considering usc given and calculating min non using pinv of system matrix

%Using pinv
xmn = pinv(A, 0.0001)*usc;

%Using SVD (min norm calculation)
xmn1 = lsqminnorm(A,usc);

%% Q18 Reshaping the solution in 2D

%Found si -> from SVD and given usc
xF = reshape(xmn, [length(x_vec), length(y_vec)]);
figure(5);
imagesc(x_vec, y_vec, abs(xF));
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
colorbar;
axis equal tight;
title(['Reconstructed constrast function of the object, M=', num2str(M)]);
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');

%% Q19 Repeating step 18 with different number of recivers and making images

%% Q20 Add noise to usc using rand and investigating influence of noise on 
%reconstruction of results

%Add noise to the usc
noiseAmp = 0.01;
uscn = usc + (((usc).*noiseAmp).*rand(length(usc), 1));

%Using pinv -> is it min norm?
xmnnoise = pinv(A)*uscn;

%Reconstruction with noise
xFn = reshape(xmnnoise, [length(x_vec), length(y_vec)]);
figure(7);
imagesc(x_vec, y_vec, abs(xFn));
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
colorbar;
axis equal tight;
title('Reconstructed constrast function of the object with noise');
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');
