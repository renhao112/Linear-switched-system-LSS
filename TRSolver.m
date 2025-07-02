function result = TRSolver(A,B,T0,TT,dt)
%the ode solver based on trapezoidal rule

nt = (TT-T0)/dt;
dn = nt/20; n = size(A,1);
result = zeros(n,20);
E = speye(n);
tempu=zeros(n,1);
[L,U] = lu(E+1/2*dt*A);
j=1;
for i=1:nt
    u=U\( L\((E-1/2*dt*A)*tempu + 1/2*dt*B*(ipt(i*dt+T0)+ipt(i*dt-dt+T0))));
    tempu = u;
    if mod(i,dn)==0
        result(:,j) = u;
        j=j+1;
    end
end


end

