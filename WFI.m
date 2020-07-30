%Final project WFT
clear;
close all;

%% Q1 and Q2 Learning about Hankel function and deriving usc
%Hankel is a third order Bessel Function given by hankle command in
%Matlab. Second question is done by substituting. 

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
% dom = [[0;0], [0;lam], [lam;lam], [lam;0]];
figure(2);
line(domX, domY, 'LineWidth', 1.5); text(lam+0.1,lam,'(\lambda,\lambda)')
hold on;
scatter(src(1), src(2), 'filled'); text(lam/2+0.1,10*lam,'  \rho_s');
% set(gca, 'YDir','reverse');
set(gca,'XAxisLocation','top','YAxisLocation','left','YDir','reverse');
grid on;
title('Figure 1: Configuration of Source and Domain in Q3');
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');

%% Q5 Defining grid and the number of grid points
step = lam/20;
x_vec = 0:step:lam;
y_vec = 0:step:lam;
[x, y] = meshgrid(x_vec, y_vec);

%Numer of points
N = length(x).*length(y);

%% Q6 Calculation of incident field on the grid
uincInit = calcUinc(x, y, src, kb);
createImage(x_vec, y_vec, uincInit);

%% Q7 Variation in background field due to change in the location of the src and kb
srcFar = [lam/2, 15*lam];
srcNear = [lam/2, 5*lam];
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
% I chose circle!! So, x^2 + y^2 = 0; Disc shape to be frank!!

% object - square 
%x_obj = lam/4:step:lam/2+lam/4;
%y_obj = lam/4:step:lam/2+lam/4;

%[x_obj,y_obj] = meshgrid(x_obj,y_obj);

%sq_eq = abs(x_obj) + abs(y_obj) - 1;

% [rho_x,rho_y] = meshgrid(g_x,g_y);

%defining krho
%k_rho = ones(size(rho_x)).*kb;
%k_rho(sq_eq <=10,sq_eq <=10) = (1.1.*kb); 

center = [lam/2, lam/2];
%center = [0, 0];
expression = (x-center(1)).^2 + (y-center(2)).^2-(lam/8).^2;
% %To make si >= 0
k_rho = ones(size(x)).*kb;
%k_rho(sq_eq <=10,sq_eq <=10) = (1.2.*kb); 
k_rho(expression <=0) = (1.2).*kb;
si_rho = (k_rho./kb).^2 - 1;

%% Q9 Making an image of the contrast function

figure(3);
imagesc(x_vec, y_vec, si_rho);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
% set(gca, 'YDir','reverse');
% view(2);
colorbar;
%caxis([min max]);
axis equal tight;
title('Contrast function of the object of choice');
xlabel('x');
ylabel('y');

%% Q10 Defining the receiver domain L = 3*lam from (-lam, 1.5lam) to (2lam, 1.5lam)
%Adding the Rec domain to figure 1

drecX = [-lam 2*lam];
drecY = [1.5*lam 1.5*lam];

figure(2);
line(drecX, drecY, 'LineWidth', 1.5);

%% Q11 Introducing grid on Drec M > 0

M = 16;
DrecX = -lam:3*lam/M:2*lam;
DrecY = 1.5*lam.*ones(size(DrecX));

figure(2);
line(DrecX, DrecY, 'LineWidth', 1.5);

%% Q12 Discretize data equation (appox. integral via sum)

%Written down in the report

%% Q13 Reshape 2D array of contrast function; study reshape function

%Reshaped contrast function
siRs = reshape(si_rho, [N, 1]);

%% Q14 Building System Matrix A from discretized data equation

%Constant in discretized equation, both x and y step are same, given in the
%question step = lam/20

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

const = -((kb^2)/16)*delX*delY*nX;

%Calculation of G/Hg, depends on 
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
title('Normalized Singular Values of System Matrix');

%s > 0 -> eigen values are non-negative -> solutions available? Read more
%on this.. 

%% Q16 Computing scattered field

usc = A*siRs;

%% Q17 Considering usc given and calculating min non using pinv of system matrix

%Using pinv -> is it min norm?
xmn = pinv(A)*usc;

%Using SVD
xmn1 = lsqminnorm(A,usc);

%% Q18 Reshaping the solution in 2D

%Found si -> from SVD and given usc
xF = reshape(xmn, [length(x_vec), length(y_vec)]);
figure(5);
imagesc(x_vec, y_vec, abs(xF));
set(gca, 'YDir','reverse');
% view(2);
colorbar;
%caxis([min max]);
axis equal tight;
title('Reconstructed constrast function of the object');
xlabel('x');
ylabel('y');

%% Q19 Repeating step 18 with different number of recivers and making images

%% Q20 Add noise to usc using rand and investigating influence of noise on 
%reconstruction of results

%Add noise to the usc
noiseAmp = 0.01;
uscn = usc + (max(max(abs(usc)))*noiseAmp).*rand(length(usc), 1);

%Using pinv -> is it min norm?
xmnnoise = pinv(A)*uscn;

%Reconstruction with noise
xF = reshape(xmnnoise, [length(x_vec), length(y_vec)]);
figure(6);
plot(abs(uscn), 'LineWidth', 1.5, 'DisplayName', 'Noise'); hold on;
plot(abs(usc), 'LineWidth', 1.5, 'DisplayName', 'Without Noise');
legend show;
figure(7);
imagesc(x_vec, y_vec, abs(xF));
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
% view(2);
colorbar;
% caxis([0 .1]);
axis equal tight;
title('Reconstructed constrast function of the object with noise');
xlabel('x');
ylabel('y');