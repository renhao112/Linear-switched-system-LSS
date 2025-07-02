function result = EBK_s(A,B,t0,tt,u0,toler,m)
%Exponential Block Krylov method
%for nonzero initial condition model
% y'=-E^{-1}Ay+E^{-1}(Bv(t)-Au0), y(0)=0, original system;
% u'= -Hu+Ev(t), u(0)=0, reduced system;
% A is stable;
% block SAI Krylov approximation without restart;

n = size (B,1); p = size(B,2)+1;
V = zeros(n  ,(m+1)*p);
H = zeros((m+1)*p,m*p);

t=tt-t0;
gamma = t/90;
E = speye(size(A));
I_gammaA = E+gamma*A;
[L,U,P,Q] = lu(I_gammaA);

%[Le,Ue,Pe,Qe] = lu(E);
B=[-A*u0,B];
%B = Qe*( Ue\( Le\(Pe*B) ) );
[Qb,Rb]=qr(B,'econ');
V(:,1:p) = Qb;

for j=1:m
    w = Q*( U\( L\(P*E*V(:,((j-1)*p+1):j*p)) ) );
    for i=1:j
        H(((i-1)*p+1):i*p,((j-1)*p+1):j*p) = V(:,((i-1)*p+1):i*p)'*w;
	    w = w - V(:,((i-1)*p+1):i*p)*H(((i-1)*p+1):i*p,((j-1)*p+1):j*p);
    end
    
    e1 = zeros(j*p,p); e1(1:p,1:p) = eye(p);
    ej = zeros(j*p,p); ej(((j-1)*p+1):j*p,1:p) = eye(p);
    invH = inv(H(1:j*p,1:j*p));
    Hjj = ( invH-eye(j*p) )/gamma;

    %reduced model
    odefun = @(t,x) -Hjj*x+e1*Rb*[1;ipt(t)];
    options = odeset('AbsTol',10^(-12),'RelTol',10^(-8)) ;
    [~,u] = ode15s(odefun,t0:t/20:tt,zeros(j*p,1),options);
    u = u';

    %residual part
    C   = I_gammaA*w;
    beta_j=zeros(n,10);
    for q=1:20
        beta_j(:,q) = C/gamma * (ej'*(invH*u(:,q+1)));
    end
    resnorm = max(vecnorm(beta_j));
    fprintf('j = %d, resnorm = %.2e\n',j,resnorm);
    if resnorm<=toler
       break
    elseif j==m
       disp('warning: no convergence within m steps');
    end

    [Q1,R1] = qr(w,'econ');
    V(:,(j*p+1):(j+1)*p) = Q1;
    H((j*p+1):(j+1)*p,((j-1)*p+1):j*p)=R1;
end

% tt = 0:100:t; nt = t/100+1;
% result = zeros(n,nt);
% for k=1:nt
%     result(:,k) = V(:,1:j)*expm(-tt(k)*Hjj)*e1*beta;
% end

result = V(:,1:j*p)*u(:,2:21)+u0*ones(1,20);


end
