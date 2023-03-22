clear variables; close all;

%Charge of the electron
q_e = 1.602e-19;

%Weight of the electron
m_e = 9.109e-31;

%Planc constant
h = 1.054572e-34;

%The value of h^2/2m
h_ = -h^2/(2.0*m_e);

%Length of the observed region
TotalLength = 5e-8;
%Resolution of the observed region
N = 400;

%Height of the potential wall confining the electron
Wallenergy = 100*q_e;
WallWidth = 10;

%Place the coordinate system
[X,Y] = meshgrid(linspace(0, TotalLength, N));
%Resolution
dx = TotalLength/N;
dy= TotalLength/N;

%a)

%Potential: V = V0 eV
V(1:N, 1:N) = 0.0;
V(1:WallWidth, :) = Wallenergy;
V(:, N-WallWidth:N) = Wallenergy;
V(:, 1:WallWidth) = Wallenergy;
V(N-WallWidth:N, :) = Wallenergy;

H = zeros(N,N);
for i = 2 : N-1
    H(i,i) = -2.0*h_/dx^2 + V(i, i);
    H(i,i-1) = h_/dx^2;
    H(i,i+1) = h_/dx^2;
end

%The sides of the matrix:
%Left side:
H(1,1) = -2.0*h_/dx^2 + V(1, 1);
H(1,2) = h_/dx^2;

%Right side:
H(N,N) = -2.0*h_/dx^2 + V(N, N);
H(N,N-1) =h_/dx^2;
[Psi,Emx] = eig(H);

%Makes a one-dimensional list from the energy matrix
E = zeros(1, N);

for i = 1:N
    E(i) = Emx(i,i);
end

%Finds the first n minimum-energy eigenvalue
n = 4;
E_minValue = zeros(n);

%Finds the minimum-energy eigenvalue
[MinE,MinEindex] = min(E);
MinE;
Psi_n(:,1) = Psi(:,MinEindex);
E_minValue(1) = MinE;

for i=2:n
    E(MinEindex) = E(MinEindex)+Wallenergy;
    [MinE,MinEindex] = min(E);
    Psi_n(:,i) = Psi(:,MinEindex);
    E_minValue(i) = MinE;
end

Hy = zeros(N,N);
for i = 2 : N-1
    Hy(i,i) = -2.0*h_/dy^2 + V(i, i);
    Hy(i,i-1) = h_/dy^2;
    Hy(i,i+1) = h_/dy^2;
end

%The sides of the matrix:
%Left side:
Hy(1,1) = -2.0*h_/dy^2 + V(1, 1);
Hy(1,2) = h_/dy^2;

%Right side:
Hy(N,N) = -2.0*h_/dy^2 + V(N, N);
Hy(N,N-1) =h_/dy^2;

[Psiy,Emxy] = eig(H);

%Makes a one-dimensional list from the energy matrix
Ey = zeros(1, N);

for i = 1:N
    Ey(i) = Emxy(i,i);
end

%Finds the first n minimum-energy eigenvalue
Ey_minValue = zeros(n);

%Finds the minimum-energy eigenvalue
[MinEy,MinEindexy] = min(Ey);
MinEy;
Psi_ny(:,1) = Psiy(:,MinEindexy);
Ey_minValue(1) = MinEy;

for i=2:n
    Ey(MinEindexy) = Ey(MinEindexy)+Wallenergy;
    [MinEy,MinEindexy] = min(Ey);
    Psi_ny(:,i) = Psiy(:,MinEindexy);
    Ey_minValue(i) = MinEy;
end

figure(11);
subplot(2, 1, 1);
mesh(X, Y ,abs(Psi_n(:,1).*transpose(Psi_ny(:,1)).^2));
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{11}|^2$', 'Interpreter',"latex");


subplot(2, 1, 2);
mesh(X, Y ,Psi_n(:,1).*transpose(Psi_ny(:,1)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{11}$', 'Interpreter',"latex");
E11 = (E_minValue(1, 1)+Ey_minValue(1, 1))/q_e

figure(12);
subplot(2, 2, 1);
mesh(X, Y ,abs(Psi_n(:,1).*transpose(Psi_ny(:,2))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{12}|^2$', 'Interpreter',"latex");

subplot(2, 2, 2);
mesh(X, Y ,Psi_n(:,1).*transpose(Psi_ny(:,2)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{12}$', 'Interpreter',"latex");
E12=(E_minValue(1, 1)+Ey_minValue(2, 1))/q_e


subplot(2, 2, 3);
mesh(X, Y ,abs(Psi_n(:,2).*transpose(Psi_ny(:,1))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{21}|^2$', 'Interpreter',"latex");

subplot(2, 2, 4);
mesh(X, Y ,Psi_n(:,2).*transpose(Psi_ny(:,1)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{21}$', 'Interpreter',"latex");
E21=(E_minValue(2, 1)+Ey_minValue(1, 1))/q_e

figure(22);
subplot(2 ,1, 1);
mesh(X, Y ,abs(Psi_n(:,2).*transpose(Psi_ny(:,2))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{22}|^2$', 'Interpreter',"latex");

subplot(2 ,1, 2);
mesh(X, Y ,Psi_n(:,2).*transpose(Psi_ny(:,2)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{22}$', 'Interpreter',"latex");
E22=(E_minValue(2, 1)+Ey_minValue(2, 1))/q_e

figure(13);
subplot(2, 2, 1);
mesh(X, Y ,abs(Psi_n(:,1).*transpose(Psi_ny(:,3))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{13}|^2$', 'Interpreter',"latex");

subplot(2, 2, 2);
mesh(X, Y ,Psi_n(:,1).*transpose(Psi_ny(:,3)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{13}$', 'Interpreter',"latex");
E13=(E_minValue(1, 1)+Ey_minValue(3, 1))/q_e

subplot(2, 2, 3);
mesh(X, Y ,abs(Psi_n(:,3).*transpose(Psi_ny(:,1))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{31}|^2$', 'Interpreter',"latex");

subplot(2, 2, 4);
mesh(X, Y ,Psi_n(:,3).*transpose(Psi_ny(:,1)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{31}$', 'Interpreter',"latex");
E31=(E_minValue(3, 1)+Ey_minValue(1, 1))/q_e

figure(15);
subplot(2, 2, 1);
mesh(X, Y ,abs(Psi_n(:,2).*transpose(Psi_ny(:,3))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{23}|^2$', 'Interpreter',"latex");


subplot(2, 2, 2);
mesh(X, Y ,Psi_n(:,2).*transpose(Psi_ny(:,3)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{23}$', 'Interpreter',"latex");
E23=(E_minValue(2, 1)+Ey_minValue(3, 1))/q_e

subplot(2, 2, 3);
mesh(X, Y ,abs(Psi_n(:,3).*transpose(Psi_ny(:,2))).^2);
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{32}|^2$', 'Interpreter',"latex");

subplot(2, 2, 4);
mesh(X, Y ,Psi_n(:,3).*transpose(Psi_ny(:,2)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{32}$', 'Interpreter',"latex");
E32=(E_minValue(3, 1)+Ey_minValue(2, 1))/q_e

figure(34);
subplot(2, 1, 1);
mesh(X, Y ,abs(Psi_n(:,3).*transpose(Psi_ny(:,3)).^2));
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{33}|^2$', 'Interpreter',"latex");

subplot(2, 1, 2);
mesh(X, Y ,Psi_n(:,3).*transpose(Psi_ny(:,3)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{33}$', 'Interpreter',"latex");
E33=(E_minValue(3, 1)+Ey_minValue(3, 1))/q_e

figure(44);
subplot(2, 1, 1);
mesh(X, Y ,abs(Psi_n(:,4).*transpose(Psi_ny(:,4)).^2));
xlabel('X');
ylabel('Y');
zlabel('$|\Psi_{44}|^2$', 'Interpreter',"latex");

subplot(2, 1, 2);
mesh(X, Y ,Psi_n(:,4).*transpose(Psi_ny(:,4)));
xlabel('X');
ylabel('Y');
zlabel('$\Psi_{44}$', 'Interpreter',"latex");
E44=(E_minValue(4, 1)+Ey_minValue(4, 1))/q_e

figure;
title('Potential');
plot3(X,Y, V/q_e,'r--');
legend('Potential (eV)');
hold off;