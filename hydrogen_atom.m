clear variables;

n=16;                   % number of solution
%Charge of the electron
q_e = -1.602e-19;
%Weight of the electron
m_e = 9.109e-31;
%Planc constant
h = 1.054572e-34;
%The value of h^2/2m
h_ = -h^2/(2.0*m_e);
%Vacuum permittivity
eps0 = 8.854187e-12;
%Length of the observed region
TotalLength = 3*5.29e-11;
%Resolution of the observed region
N = 120;
k_e=1/(4*pi*eps0);

Wallenergy = -300*q_e;

x=linspace(-TotalLength/2,TotalLength/2,N);
y=linspace(-TotalLength/2,TotalLength/2,N);

%Place the coordinate system in the middle (0 is in the middle)
[X,Y] = meshgrid(x, y);
%Resolution
dx = TotalLength/N;
dy= TotalLength/N;


R=5.29e-11;      % r_x
x0=0; y0=0;      % center

r0=5.29e-11;

idx=((X-x0)/R).^2 + ((Y-y0)/R).^2 < 1;

Vb=-300*q_e;                 % Potential barrier

V=zeros(N, N);
Vp=zeros(N, N);

for i=1:N
    for j=1:N
        r_=sqrt(x(i)^2+y(j)^2);
            if r_==0
                Va=(-k_e*q_e^2)/(1e-12);
            else
                Va=-k_e*q_e^2/r_;
            end
        Vp(i, j)= Va;
    end
end

V=(idx)*0 + (1-idx)*Vb;

V=V+Vp;

Dx = zeros(N, N);

for i = 2 : N-1
    Dx(i,i) = -2.0*h_/dx^2;
    Dx(i,i-1) = 1*h_/dx^2;
    Dx(i,i+1) = 1*h_/dx^2;
end

a1=1;
a2=0.95;
%The sides of the matrix:
%Left side:
Dx(1,1) = -2*h_/dx^2;
Dx(1,2) = 1*h_/dx^2;

%Right side:
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

%Right side:
Dy(N, N)= -2*h_/dy^2;
Dy(N, N-1)= 1*h_/dy^2;

Dy=sparse(Dy);

%kronsum
T =-0.5*((kron(Dy, speye(N))) + kron(speye(N), Dx));

T = sparse(T);

U=diag(reshape(V, N^2, 1));
U = sparse(U);

H=T+U;
[Psi,Emx] = eigs(H, N);

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

Energy= transpose(E_minValue(1:6, 1))./(-q_e) %meV

figure(4);
for i=1:6
    subplot(2, 3, i);
    pcolor(X, Y, abs(Psi_n(:,:,i)).^2);
    xlabel('X')
    ylabel('Y')
    %axis equal
    shading flat
    colormap(jet)
end

Energy=transpose(E_minValue(7:12, 1))./(-q_e) %eV

figure(5);
for i=7:12
    subplot(2, 3, i-6);
    pcolor(X, Y, abs(Psi_n(:,:,i)).^2);
    xlabel('X')
    ylabel('Y')
    %axis equal
    shading flat
    colormap(jet)
end

Energy=transpose(E_minValue(13:16, 1))./(-q_e) %eV

figure(6);
for i=13:16
    subplot(2, 2, i-12);
    pcolor(X, Y, abs(Psi_n(:,:,i)).^2);
    xlabel('X')
    ylabel('Y')
    %axis equal
    shading flat
    colormap(jet)
end