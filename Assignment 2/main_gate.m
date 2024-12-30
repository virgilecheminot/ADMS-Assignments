clear; close all; clc

%% Question 1
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('ADMS_FEM_05');

dis_stru(posit,l,gamma,xy,pr,idb,ndof);

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Question 2
[MFF, MFC, MCC] = freefree(M,ndof);
[KFF, KFC, KCC] = freefree(K,ndof);

[x0, omega_squared] = eig(MFF\KFF);
omega = diag(sqrt(omega_squared));
[omega, ind] = sort(omega);
modes = x0(:,ind);

sf = 2;
figure
t = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ii=1:4
    nexttile
    diseg2(modes(:,ii), sf, incid, l, gamma, posit, idb, xy, 0.12);
    title(sprintf('Mode %d - Frequency: %.2f Hz', ii, omega(ii)/(2*pi)))
    legend('Location', 'South')
end

%% Question 3
B = [0.02 0.03]';
A = zeros(2,2);
for ii=1:2
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
ab = A\B;
C = ab(1)*M + ab(2)*K;
[CFF, CFC, CCC] = freefree(C, ndof);

F0 = zeros(ndof,1);
F0(idb(11,1)) = 1;
om = (0:0.01:10)*2*pi;

X = zeros(ndof, length(om));
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

%% Question 4
ii = 1:2;
Phi = modes(:,ii);
Mm = Phi'*MFF*Phi;
Km = Phi'*KFF*Phi;
Cm = Phi'*CFF*Phi;
Fm = Phi'*F0;

Xmod = zeros(2, length(om));
for jj=1:length(om)
    A = -om(jj)^2*Mm + 1i*om(jj)*Cm + Km;
    Xmod(:,jj) = A\Fm;
end
Xm = Phi*Xmod;

figure();
subplot(2,1,1);
semilogy(om/(2*pi), abs(X(idb(9,1),:)), 'LineWidth', 1.5);
hold on;
semilogy(om/(2*pi), abs(Xm(idb(9,1),:)), 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('Amplitude [m/N]');
title('Frequency response function');
grid on;
legend('Original', 'Modal superposition', 'Location', 'southwest');
subplot(2,1,2);
plot(om/(2*pi), angle(X(idb(9,1),:)), 'LineWidth', 1.5);
hold on;
plot(om/(2*pi), angle(Xm(idb(9,1),:)), 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('Phase [rad]');
grid on;

%% Question 6
p0 = -1e2; % N/m
p0G = [0 p0]';
FG = zeros(3*nnod,1);

for ii = 10:15
    Lk = l(ii);
    g = gamma(ii);
    p0L = [cos(g) sin(g);
          -sin(g) cos(g)]*p0G;
    FkL = [Lk/2 ;  0   ; 0 ; Lk/2 ;  0   ; 0] * p0L(1) + ...
          [ 0   ; Lk/2 ; 0 ;  0   ; Lk/2 ; 0] * p0L(2);
    FkG = [cos(g) -sin(g) 0 0       0      0;
           sin(g)  cos(g) 0 0       0      0;
           0       0      1 0       0      0;
           0       0      0 cos(g) -sin(g) 0;
           0       0      0 sin(g)  cos(g) 0;
           0       0      0 0       0      1] * FkL;
    Ek = zeros(6,3*nnod);
    for jj = 1:6
        Ek(jj,incid(ii,jj)) = 1;
    end
    FG = FG + Ek'*FkG;
end

xF = KFF\FG(1:ndof);
figure();
diseg2(xF, 20, incid, l, gamma, posit, idb, xy, 0.05);
txt = sprintf('\\Delta y_1 = %.4f m', xF(idb(1,2)));
text(0,5,txt,'HorizontalAlignment','center', 'FontSize', 12);
title('Static deflection')
legend('Location', 'South')