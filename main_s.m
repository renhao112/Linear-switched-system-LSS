%to get the CPU time for serial computation

clear; clc;

%% initialize the system
Nx = 400; n = Nx*Nx;
[As,Bs] = generate_lss(Nx);
Nt = 100; load('refsol4.mat');
%% serial block Krylov subspace solution
cptime_s4 = zeros(3,1);
erors4 = zeros(3,Nt);
for k = 1:3
    toler = 10^(-(k+2));
    ssol = zeros(n,Nt);
    tic
    ssol(:,1:20) = EBK(As{1},Bs{1},0,2,toler,50);
    for i = 2:5
        ssol(:,(1:20)+(i-1)*20) = EBK_s(As{i},Bs{i},2*(i-1),2*i,ssol(:,(i-1)*20),toler,50);
    end
    cptime_s4(k) = toc;
    erors4(k,:) = vecnorm(ssol-refsol(:,2:Nt+1))./vecnorm(refsol(:,2:Nt+1));
end

%% serial Tr solution
dt = 0.5*10^(-2); 
trsol = zeros(n,Nt+1);
tic
for i = 1:5
    trsol(:,(i-1)*20+(1:20)+1)=TR_s(As{i},Bs{i},2*(i-1),2*i,trsol(:,(i-1)*20+1),dt);
end
cptime_trs4 = toc;
erortrs4 = vecnorm(trsol(:,2:Nt+1)-refsol(:,2:Nt+1))./vecnorm(refsol(:,2:Nt+1));