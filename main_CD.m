%artificial linear switched system from a convection diffusion equation
% u_t - \nu\Delta u + b\cdot\grad u = f, x in [0,1]x[0,1]
%homogeneous boundary condition, zero initial value condition
%spatial discretization: central differential scheme, gives
% x_t + As x = Bs u, switching piont t=2,4,6,8

clear; clc;

%% initialize the system
Nx = 300; n = Nx*Nx;
[As,Bs] = generate_lss(Nx);

%% generate the reference solution
refsol = zeros(n,20*5+1);
%refsol(:,1) = u0;
refsol(:,(1:20)+1) = EBK(As{1},Bs{1},0,2,10^(-10),70);
for i = 2:5
    refsol(:,(1:20)+(i-1)*20+1) = EBK_s(As{i},Bs{i},2*(i-1),2*i,refsol(:,(i-1)*20+1),10^(-10),70);
end

%% paraexp krylov
toler = 10^(-5); Nt = 100;
ftoler = toler/10;
parasol = zeros(n,Nt);
cptime_35 = zeros(5,1);

%EBK solver
% tempsol = zeros(n,Nt);
% tic
% tempsol(:,1:20)=EBK_s(As{1},Bs{1},0,2,u0,toler,50);
% %SAI solver
% for j=2:5
%    tempsol(:,(j-1)*20+(1:20))=SAI_appro(As{j},tempsol(:,(j-1)*20),2,ftoler,30);
% end
% cptime_4(1) = toc;
% parasol=parasol+tempsol;

for i=1:5
    tempsol = zeros(n,Nt);
    %EBK solver
    tic
    tempsol(:,(i-1)*20+(1:20))=EBK(As{i},Bs{i},2*(i-1),2*i,toler,50);
    %SAI solver
    for j=i+1:5
       tempsol(:,(j-1)*20+(1:20))=SAI_appro(As{j},tempsol(:,(j-1)*20),2,ftoler,40);
    end
    cptime_35(i) = toc;
    parasol=parasol+tempsol;
    %eror(i,:) = vecnorm(parasol-refsol(:,2:401));
end
%errot_5 = vecnorm(C*(parasol-refsol(:,2:401)))./vecnorm(C*refsol(:,2:401));
erorst_35 = vecnorm(parasol-refsol(:,2:Nt+1))./vecnorm(refsol(:,2:Nt+1));


%% paraexp Tr
dt = 0.5*10^(-2); Nt = 100;
%ftoler = toler/10;
trsol = zeros(n,Nt);
cptime_tr3 = zeros(5,1);

%EBK solver

for i=1:5
    tempsol = zeros(n,Nt);
    %Tr solver
    tic
    tempsol(:,(i-1)*20+(1:20))=TRSolver(As{i},Bs{i},2*(i-1),2*i,dt);
    %SAI solver
    for j=i+1:5
       tempsol(:,(j-1)*20+(1:20))=SAI_appro(As{j},tempsol(:,(j-1)*20),2,10^(-5),40);
    end
    cptime_tr3(i) = toc;
    trsol=trsol+tempsol;
    %eror(i,:) = vecnorm(parasol-refsol(:,2:401));
end

erortr3 = vecnorm(trsol-refsol(:,2:Nt+1))./vecnorm(refsol(:,2:Nt+1));
%%
xspan = (1:100)/10;
semilogy(xspan,erorst_14,'Color',"#D95319",'LineStyle',':','LineWidth',1.2);
hold on
semilogy(xspan,erortr4,'Color',"#77AC30",'LineStyle','--','LineWidth',1.2);
hold on
%% plot
xspan = (1:100)/10;
semilogy(xspan,erorst_33,'Color',"#D95319",'LineStyle',':','LineWidth',1.2);
hold on
semilogy(xspan,erorst_34,'Color',"#77AC30",'LineStyle','--','LineWidth',1.2);
hold on
semilogy(xspan,erorst_35,'Color',"#EDB120",'LineStyle','-','LineWidth',1.2);
hold on
semilogy(xspan,erortr3,'Color',"#0072BD",'LineStyle','--','LineWidth',1.2);