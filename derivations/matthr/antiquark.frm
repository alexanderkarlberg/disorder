* Tree-level |M|^2 structure of  l(k) qbar(p1) -> l(kp) qbar(p2) g(p3),
* for each lepton helicity l and chirality h of the quark field
* (6 = right-handed, 7 = left-handed); the projector acts on v(p2).
* Fermion line: vbar(p1) [ g^mu (-p2-p3) g^rho /(p2+p3)^2
*                        + g^rho (p3-p1) g^mu /(p1-p3)^2 ] v(p2)
Vectors k,kp,p1,p2,p3;
Indices mu,nu,rho;
Symbols is13,is23,S13,S23;
#do l = 6,7
#do h = 6,7
Local X`l'`h' = S13^2*S23^2*
   (1/2)*g_(1,kp)*g_(1,mu)*g`l'_(1)*g_(1,k)*g_(1,nu)
 * (1/2)*(-1)*
     g_(2,p1)*( g_(2,mu)*(-g_(2,p2)-g_(2,p3))*g_(2,rho)*is23
              + g_(2,rho)*(g_(2,p3)-g_(2,p1))*g_(2,mu)*(-is13) )
     *g`h'_(2)*g_(2,p2)*
            ( g_(2,rho)*(-g_(2,p2)-g_(2,p3))*g_(2,nu)*is23
            + g_(2,nu)*(g_(2,p3)-g_(2,p1))*g_(2,rho)*(-is13) );
#enddo
#enddo
trace4,1;
trace4,2;
.sort
repeat id is13*S13 = 1;
repeat id is23*S23 = 1;
id S13 = 2*p1.p3;
id S23 = 2*p2.p3;
id p3 = k + p1 - kp - p2;
id p1.p1 = 0; id p2.p2 = 0; id k.k = 0; id kp.kp = 0;
id p1.p2 = k.p1 - k.kp - k.p2 - kp.p1 + kp.p2;
.sort
Format 255;
Format nospaces;
Print;
.end
