* Tree-level |M|^2 structure of  l(k) g(p1) -> l(kp) q(p2) qbar(p3),
* for each lepton helicity l and quark-line chirality h
* (6 = right-handed, 7 = left-handed), labels as in DISENT.
* Couplings, colour, 1/q^4 stripped, gluon polarisation sum -g.
* X(l,h) = A(l,h) * s12^2 * s13^2 (s12 = 2 p1.p2, s13 = 2 p1.p3).
Vectors k,kp,p1,p2,p3;
Indices mu,nu,rho;
Symbols is12,is13,S12,S13;
#do l = 6,7
#do h = 6,7
Local X`l'`h' = S12^2*S13^2*
   (1/2)*g_(1,kp)*g_(1,mu)*g`l'_(1)*g_(1,k)*g_(1,nu)
 * (1/2)*(-1)*
     g_(2,p2)*( g_(2,rho)*(g_(2,p2)-g_(2,p1))*g_(2,mu)*(-is12)
              + g_(2,mu)*(g_(2,p1)-g_(2,p3))*g_(2,rho)*(-is13) )
     *g`h'_(2)*g_(2,p3)*
            ( g_(2,nu)*(g_(2,p2)-g_(2,p1))*g_(2,rho)*(-is12)
            + g_(2,rho)*(g_(2,p1)-g_(2,p3))*g_(2,nu)*(-is13) );
#enddo
#enddo
trace4,1;
trace4,2;
.sort
repeat id is12*S12 = 1;
repeat id is13*S13 = 1;
id S12 = 2*p1.p2;
id S13 = 2*p1.p3;
id p3 = k + p1 - kp - p2;
id p1.p1 = 0; id p2.p2 = 0; id k.k = 0; id kp.kp = 0;
id p1.p2 = k.p1 - k.kp - k.p2 - kp.p1 + kp.p2;
.sort
Format 255;
Format nospaces;
Print;
.end
