from operator import mul
from functools import reduce
pi=3.141592653589793
e=2.718281828459045
ln2=0.6931471805599453
euler_gamma=0.5772156649015329

#One line integration (rectangular)
integrate=lambda f,a,b,dx=0.0001:sum(map(lambda x:f(x*dx+a)*dx,range(int(b/dx)-int(a/dx)))) if a<b else -integrate(f,b,a,dx)

#One line derivate
derivate=lambda f,x,dx=0.0001,side=0:((f(x+dx)-f(x))/dx + (f(x-dx)-f(x))/-dx)/2 if side==0 else (f(x+dx)-f(x))/dx if side>0 else (f(x-dx)-f(x))/-dx

#One line factorial
factorial=lambda n:reduce(mul,range(2,n+1),1)

#One line pgcd
pgcd=lambda a,b:abs(a) if b==0 else pgcd(abs(b),abs(a)%abs(b))

#One line ppcm
ppcm=lambda a,b:abs(a*b)//pgcd(a,b)

#One line sin DL11 approximation
sin=lambda x:(x%(2*pi)-pi)**11/39916800-(x%(2*pi)-pi)**9/362880+(x%(2*pi)-pi)**7/5040-(x%(2*pi)-pi)**5/120+(x%(2*pi)-pi)**3/6+pi-x%(2*pi)

#One line arithmetic and geometric means
aagm=lambda x,y,n=10:0.5*(x+y) if n==1 else 0.5*(aagm(x,y,n-1)+gagm(x,y,n-1))
gagm=lambda x,y,n=10:(x*y)**0.5 if n==1 else (aagm(x,y,n-1)*gagm(x,y,n-1))**0.5

#One line ln approximation
ln=lambda x,m=50:pi/(2*aagm(1,4/(x*2**m)))-m*ln2

#One line li approximation
li=lambda x,m=50:euler_gamma+ln(ln(x))+x**0.5*sum(map(lambda n:((-1)**(n-1)*ln(x)**n)/(factorial(n)*2**(n-1))*sum(map(lambda k:1/(2*k+1),range(int((n-1)/2)+1))),range(1,m+1)))

#One line gamma Lanczo approximation
gamma=lambda z:((2*pi)**0.5)*((complex(z)+6.5)**(complex(z)-0.5))*(e**(-complex(z)-6.5))*sum(map(lambda i:[0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7][i]/((complex(z)-1+i) if i!=0 else 1),range(9))) if complex(z).real>=0.5 else pi/(sin(pi*z)*gamma(1-z))

#One line zeta alternating series approximation
dk=lambda k,n:n*sum(map(lambda j:(factorial(n+j-1)*4**j)/(factorial(n-j)*factorial(2*j)),range(k,n+1)))
zeta=lambda s,n=50:1/(dk(0,n)*(1-2**(1-s)))*sum(map(lambda k:(-1)**(k-1)*dk(k,n)/(k**s),range(1,n+1)))
