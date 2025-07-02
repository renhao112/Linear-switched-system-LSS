function result= SAI_appro(A,v,t,toler,m)
%SAI approximation for exp(-sA)v, s in [0,t]

n = size (v,1);
V = zeros(n  ,m+1);
H = zeros(m+1,m);

E = speye(size(A));
gamma = t/10;
I_gammaA = E+gamma*A;
%fprintf('Computing sparse LU factorization of the SAI matrix...');
[L,U,P,Q] = lu(I_gammaA);

beta = norm(v);
V(:,1) = v/beta;
options = odeset('AbsTol',10^(-12),'RelTol',10^(-8)) ;

for j=1:m
    w = Q*( U\( L\(P*E*V(:,j)) ) );
    for i=1:j
        H(i,j) = w'*V(:,i);
	w      = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    e1 = zeros(j,1); e1(1) = 1;
    ej = zeros(j,1); ej(j) = 1;
    invH = inv(H(1:j,1:j));
    Hjj = ( invH-eye(j,j) )/gamma;
    C   = norm(I_gammaA*w);
    s   = (1/9:1/9:1)*t;
    for q=1:length(s)  %expm solver
        u         = expm(-s(q)*Hjj)*e1*beta;
        beta_j(q) = C/gamma * (ej'*(invH*u));
    end
    resnorm = norm(beta_j,'inf');
    fprintf('j = %d, resnorm = %.2e\n',j,resnorm);
    if resnorm<=toler
       break
    elseif j==m
       disp('warning: no convergence within m steps');
    end
    V(:,j+1) = w/H(j+1,j);
end

tt = 0:t/20:t; 
result = zeros(n,20);
for k=1:20
    result(:,k) = V(:,1:j)*expm(-tt(k+1)*Hjj)*e1*beta;
end


end

