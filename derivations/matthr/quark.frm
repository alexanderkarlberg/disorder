* Tree-level |M|^2 structure of  l(k) q(p1) -> l(kp) q(p2) g(p3),
* exchanged boson momentum q = k - kp = p2 + p3 - p1,
* for each lepton helicity l and quark helicity h
* (6 = right-handed projector (1+g5), 7 = left-handed (1-g5)).
* Momentum labels as in DISENT: 1 = incoming parton, 2,3 = outgoing
* partons, k = P(,6) incoming lepton, kp = P(,7) outgoing lepton.
* Couplings, colour factors, 1/q^4 are stripped:
*   A(l,h) = L_{mu nu} H^{mu nu},  gluon polarisation sum -g_{rho sigma}
* X(l,h) = A(l,h) * s13^2 * s23^2 (s13 = 2 p1.p3, s23 = 2 p2.p3) is a
* polynomial; output in terms of independent invariants with
* p3 = k + p1 - kp - p2 eliminated, and p1.p2 eliminated via p3^2 = 0.
Vectors k,kp,p1,p2,p3;
Indices mu,nu,rho;
Symbols is13,is23,S13,S23;
#do l = 6,7
#do h = 6,7
Local X`l'`h' = S13^2*S23^2*
   (1/2)*g_(1,kp)*g_(1,mu)*g`l'_(1)*g_(1,k)*g_(1,nu)
 * (1/2)*(-1)*
     g_(2,p2)*( g_(2,rho)*(g_(2,p2)+g_(2,p3))*g_(2,mu)*is23
              - g_(2,mu)*(g_(2,p1)-g_(2,p3))*g_(2,rho)*is13 )
     *g`h'_(2)*g_(2,p1)*
            ( g_(2,nu)*(g_(2,p2)+g_(2,p3))*g_(2,rho)*is23
            - g_(2,rho)*(g_(2,p1)-g_(2,p3))*g_(2,nu)*is13 );
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
