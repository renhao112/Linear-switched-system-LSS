unction [As,Bs] = generate_lss(Nx)
%artificial linear switched system from a convection diffusion equation
% u_t - \nu\Delta u + b\cdot\grad u = f, [0,1]x[0,1]
%homogeneous boundary condition, zero initial value condition
%spatial discretization: central differential scheme, gives
% x_t + A x = F{\sigma} u, switching piont t=2,4,6

dx = 1/(Nx+1);
nu = [0.06,0.04,0.02,0.08,0.1];
K1 = 1/(dx)^2*spdiags([-ones(Nx,1),2*ones(Nx,1),-ones(Nx,1)],-1:1,Nx,Nx);
B1 = 1/(2*dx)*spdiags([-ones(Nx,1),zeros(Nx,1),ones(Nx,1)],-1:1,Nx,Nx);
xx = dx:dx:Nx*dx;
D1 = spdiags((1-xx.^2)',0,Nx,Nx); D2 = spdiags(2*xx',0,Nx,Nx);

B = zeros(Nx*Nx,10);
for i = 1:10
    b = zeros(1,Nx); b(i:10:Nx) = 1; B(:,i) = reshape(kron(b,ones(Nx,1)),[],1)/Nx;
end
Bs = {B(:,[1,3,4,7]),B(:,[2,3,4,8]),B(:,[4,5,6,9]),B(:,[7,8,9,10]),B(:,[7,3,5,6])};

As = cell(5,1);
for i = 1:5
    As{i} = kron(nu(i)*K1,speye(Nx))+kron(speye(Nx),nu(i)*K1)+kron(D2,D1*B1)+kron(D1*B1,-D2); 
end

%u0 = reshape(sin(pi*xx')*sin(pi*xx),[],1);

end