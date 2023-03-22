clear variables;
n=18;                   % number of solution
%Charge of the electron
q_e = 1.602e-19;
%Weight of the electron
m_e = 9.109e-31;
%Planc constant
h = 1.054572e-34;
%The value of h^2/2m
h_ = h^2/(2.0*m_e);
%Vacuum permittivity
%Length of the observed region
TotalLength = 5e-8;
%Resolution of the observed region
N = 70;

Wallenergy = 100*q_e;
WallWidth = 10;
x=linspace(0,TotalLength,N);
y=linspace(0,TotalLength,N);

%Place the coordinate system
[X,Y] = meshgrid(x, y);
%Resolution
dx = TotalLength/N;
dy= TotalLength/N;


V(1:N, 1:N) = 0.0;
V(1:WallWidth, :) = Wallenergy;
V(:, N-WallWidth:N) = Wallenergy;
V(:, 1:WallWidth) = Wallenergy;
V(N-WallWidth:N, :) = Wallenergy;

Dx = zeros(N, N);
for i = 2 : N-1
    Dx(i,i) = -2.0*h_/dx^2;
    Dx(i,i-1) = 1*h_/dx^2;
    Dx(i,i+1) = 1*h_/dy^2;
end
a1=1;
a2=0.95;
%The sides of the matrix:
%Left side:
Dx(1,1) = -2*h_/dx^2;
Dx(1,2) = 1*h_/dx^2;
%Right side: Psi(PrSize+1)=0
Dx(N, N)= -2*h_/dx^2;
Dx(N, N-1)= 1*h_/dx^2;

Dy = zeros(N, N);
for i = 2 : N-1
    Dy(i,i) = -2.0*h_/dy^2;
    Dy(i,i-1) = 1*h_/dy^2;
    Dy(i,i+1) = 1*h_/dy^2;
end

%The sides of the matrix:
%Left side:
Dy(1,1) = -2*h_/dy^2;
Dy(1,2) = 1*h_/dy^2;
%Right side: Psi(PrSize+1)=0
Dy(N, N)= -2*h_/dy^2;
Dy(N, N-1)= 1*h_/dy^2;

%kronsum
%T=-0.5*((kron(D,eye(size(D,2))) + kron(eye(size(D,2)),D)));
T = -0.5*((a1*kron(Dy, speye(N))) + a2*kron(speye(N), Dx)) ;
U=diag(reshape(V, N^2, 1));
H=T+U;
[Psi,Emx] = eig(H);

%Makes a one-dimensional list from the energy matrix
E = zeros(1, N);
for i = 1:N
    E(i) = Emx(i,i);
end
%Finds the first n minimum-energy eigenvalue
E_minValue = zeros(n);
%Finds the minimum-energy eigenvalue
[MinE,MinEindex] = min(E);
MinE;
Psi_n(:, :,1) = reshape(transpose(Psi(:,MinEindex)), N, N);
E_minValue(1) = MinE;
for i=2:n
    E(MinEindex) = E(MinEindex)+Wallenergy;
    [MinE,MinEindex] = min(E);
    Psi_n(:,:,i) = reshape(transpose(Psi(:,MinEindex)), N, N);
    E_minValue(i) = MinE;
end

Energy=transpose(E_minValue(1:9, 1)./q_e)
figure(2);
for i=1:9
    subplot(3,3, i);
    mesh(X, Y, Psi_n(:,:,i));
    xlabel('x');
    ylabel('y');
    zlabel('$\Psi$', 'Interpreter',"latex");
end

Energy=transpose(E_minValue(10:18, 1)./q_e)

figure(3);
for i=10:18
    subplot(3,3, i-9);
    mesh(X, Y, Psi_n(:,:,i));
    xlabel('x');
    ylabel('y');
    zlabel('$\Psi$', 'Interpreter',"latex");
end

figure(9);

L=1;
x=linspace(0,L,N);
y=linspace(0,L,N);

[X,Y] = meshgrid(x, y);

%Analytical solution:

subplot(3,3,1);
PSI11= 2/((L^2)^0.5).*sin((1*pi/L).*X).*sin((1*pi/L).*Y);
mesh(X, Y, PSI11);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,2);
PSI12= 2/((L^2)^0.5).*sin((1*pi/L).*X).*sin((2*pi/L).*Y);
mesh(X, Y, PSI12);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,3);
PSI21= 2/((L^2)^0.5).*sin((2*pi/L).*X).*sin((1*pi/L).*Y);
mesh(X, Y, PSI21);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,4);
PSI13= 2/((L^2)^0.5).*sin((1*pi/L).*X).*sin((3*pi/L).*Y);
mesh(X, Y, PSI13);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,5);
PSI31= 2/((L^2)^0.5).*sin((3*pi/L).*X).*sin((1*pi/L).*Y);
mesh(X, Y, PSI31);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,6);
PSI23= 2/((L^2)^0.5).*sin((2*pi/L).*X).*sin((3*pi/L).*Y);
mesh(X, Y, PSI23);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,7);
PSI32= 2/((L^2)^0.5).*sin((3*pi/L).*X).*sin((2*pi/L).*Y);
mesh(X, Y, PSI32);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,8);
PSI33= 2/((L^2)^0.5).*sin((3*pi/L).*X).*sin((3*pi/L).*Y);
mesh(X, Y, PSI33);
zlabel('$\Psi$', 'Interpreter',"latex");

subplot(3,3,9);
PSI44= 2/((L^2)^0.5).*sin((4*pi/L).*X).*sin((4*pi/L).*Y);
mesh(X, Y, PSI44);
zlabel('$\Psi$', 'Interpreter',"latex");