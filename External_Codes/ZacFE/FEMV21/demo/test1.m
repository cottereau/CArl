%This is an example of an incident EM-wave in vacuum hitting
%a glass-plate of thickness 2a. The glassplate has the conductivity
%sigma=.02 mho/m and the relative permittivity of e_r=2.5.

if exist('w')~=1,w=30e9;end
if exist('n')~=1,n=7;end
ww=input(['w [' num2str(w/1e9) ' GHz]: '])*1e9;
if ~isempty(ww),w=ww;end
nn=input(['n (# refinements) [' num2str(n) ']: ']);
if ~isempty(nn),n=nn;end
fprintf(['\n# elements: ' num2str(2*2^n) '\n'])
a=.02;b=2*a;
c0=299792458;
mu0=pi*4e-7;
e0=1/(c0^2*mu0);

xn=refine1([-b -a a b],n);
N=size(xn);
alpha=ones(N);
eps=e0*[(abs(xn)>a)+2.5*(abs(xn)<=a)];
sigma=.02*(abs(xn)<a);
beta=mu0*[j*w*sigma-w^2*eps];
s=zeros(N);

k0=w/c0;
E0=1;
Ei=-2j*k0*E0*exp(-j*k0*xn(1));
Ef=0;
gi=-j*k0;
gf=j*k0;
bound=...
   [nan Ei gi
    nan Ef gf]
z=fem1(xn,alpha,beta,s,bound);
plot(xn,abs(z),'.-')
grid