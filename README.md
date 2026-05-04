# Dokumentáció
A projekt lényegi része a contrib mappában található a fontosabb fájlok:

- test_run.cc, test_runs.cc: Futatásra használt eszközök, minkettő determinsztikus random generált problémákat old meg. A testrun argumentumai: n(int), seed(int), Mit futasson(stringek szóközzel elválasztva) BnC HK Heu Log, ezek az opciók. A testrun tetszőleges számú pozitív egész fogad el, ezekre futatja a beállított feladatokat. Példa futások:
  - time ./contrib/test_run 30 2 BnC Heu  | Ez megoldja a 30 csúcsú 2es seed feladatot a BnC illetve a Heu algoritmussal.
  - time ./contrib/test_runs 10 20 30 40  | Ez megold négy feladatot. A forráskódban a whattodo vector állítássával lehet befolyásolni mit végez el
 
- BnC.h: A vágósíkos algoritmus

- Heldkarp.h: A Held-Karp továbbfejlesztése, jelenleg nem teljesen működőképes

- Heuristic.cc Néhány heurisztika

