function u = ipt(t)
%input function
%u1 = 1/2sin(t/20)e^{-t/500}+1/20e^{-t/500};
%u2 = u3 = 1/2, constant;
%u4, u5, u6

u=zeros(4,1);
u(1) = -sin(6*pi*t)*exp(-t/2);
u(2) = 2*cos(4*pi*t);
u(3) = abs(4*sin(2*pi*t))*exp(-t);
u(4) = -3;

end

