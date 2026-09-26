import sys, subprocess
from evalprog import values
pts = sys.argv[1]
G = values("qgg_g", pts); DA = values("d_a", pts); DB = values("d_b", pts); E = values("e", pts)
for CF, CA in ((4/3, 3.0), (0.7, 2.2)):
    out = subprocess.run(["harness/harness", "matfor", pts, str(CF), str(CA), "0.5"], capture_output=True, text=True).stdout.split("\n")
    print(f"CF={CF:.3f} CA={CA}")
    for g, da, db, e, line in zip(G, DA, DB, E, out):
        QG, GG, D1, D2, EE, EMSQ = map(float, line.split())
        gm = CF*(g["X1h77"]+g["X1h76"]) + (CF-CA/2)*(g["X2h77"]+g["X2h76"])
        print(f"  G {GG/gm:.12f}  D1/d_a {D1/(da['Dh77']+da['Dh76']):.12f}  D1/d_b {D1/(db['Dh77']+db['Dh76']):.6f}"
              f"  D2/d_b {D2/(db['Dh77']+db['Dh76']):.12f}  E {EE/(e['Eh77']+e['Eh76']):.12f}")
