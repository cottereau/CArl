x0=linspace(-1,1,10);
alpha=inline('x | 1');
beta=inline('x | 1');
source=inline('x & 0');
bd=[-1 nan nan;1 nan nan];

u=inline('exp(1)/(1-exp(2))*(exp(-x)-exp(x))');
tol=input('tol: ');
[x1,a1,b1,s1]=adapt1(x0,tol,alpha,beta,source,bd);
U1=fem1(x1,a1,b1,s1,bd).';
figure,plot(x1,U1,'.-',x1,u(x1))
figure,plot(x1,abs(u(x1)-U1));
figure,plot(diff(x1),'.-')