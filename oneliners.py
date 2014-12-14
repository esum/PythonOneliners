from math import sin
pi=3.141592653589793
e=2.718281828459045

#One line integration (rectangular)
integrate=lambda f,a,b,dx=0.0001:sum(map(lambda x:f(x*dx+a)*dx,range(0,int(b/dx)-int(a/dx))))

#One line gamma Lanczo approximation
gamma=lambda z:((2*pi)**0.5)*((complex(z)+6.5)**(complex(z)-0.5))*(e**(-complex(z)-6.5))*sum(map(lambda i:[0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7][i]/((complex(z)-1+i) if i!=0 else 1),range(0,9))) if complex(z).real>=0.5 else pi/(sin(pi*z)*gamma(1-z))
