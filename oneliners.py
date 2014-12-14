from operator import mul
from functools import reduce
pi=3.141592653589793
e=2.718281828459045

#One line integration (rectangular)
integrate=lambda f,a,b,dx=0.0001:sum(map(lambda x:f(x*dx+a)*dx,range(0,int(b/dx)-int(a/dx))))

#One line factorial
factorial=lambda n:reduce(mul,range(2,n+1),1)

#One line sin DL11
sin=lambda x:(x%(2*pi)-pi)**11/39916800-(x%(2*pi)-pi)**9/362880+(x%(2*pi)-pi)**7/5040-(x%(2*pi)-pi)**5/120+(x%(2*pi)-pi)**3/6+pi-x%(2*pi)

#One line gamma Lanczo approximation
gamma=lambda z:((2*pi)**0.5)*((complex(z)+6.5)**(complex(z)-0.5))*(e**(-complex(z)-6.5))*sum(map(lambda i:[0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7][i]/((complex(z)-1+i) if i!=0 else 1),range(0,9))) if complex(z).real>=0.5 else pi/(sin(pi*z)*gamma(1-z))

#One line zeta alternating series approximation
dk=lambda k,n:n*sum(map(lambda j:(factorial(n+j-1)*4**j)/(factorial(n-j)*factorial(2*j)),range(k,n+1)))
zeta=lambda s,n=50:1/(dk(0,n)*(1-2**(1-s)))*sum(map(lambda k:(-1)**(k-1)*dk(k,n)/(k**s),range(1,n+1)))
