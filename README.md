# Thesis
## Todo
- #### Működés
  - A solver osztály templatté alakítása
  - Label indexelésének ellenőrzése, helyetesítése
  - Ha nincs megoldás viselkedés
  - Nem csak 0 kezdetű utak?
- #### Optimalizálás
  - AP contractciós becslés
  - Az éleket egy taskon (vagy az egész alatt) iterálásban csökkenő sorrendben adjuk hozzá.
    - ehez kupacok?
      - Operátorok a task struktúrába
  - Upper bound legyen élköltség összeg (eközben gyüjthetünk más infot is) a std::limits szám helyett, vagy fusson le először egy közelítő verzió ami a korlátot generál, közben néhány maszkra alsó korlátot generál, potenciálisan iterált verzió?
- #### Debug
  - Warningok átnézése
  - Esetleg néhány bound függvényhez Inline kulcsszó hozzáadása?
- #### Tesztelés
  - A Rand.cc osztályba szervezése
  - Nagy teszt program írása
<br>

## Fájlok
- RndGraph.cc - Ez a random gráfok generálását végző script. tartalmaz random példa generáló eszközöket amiket az elöző félév során teszteléshez írtam esetleg jövőben hasznosak lehetnek.
- smallTSPsolver.cc - A kis (egyenlőre csak zárt, 0-ból induló) TSP feladatok megoldására szolgáló osztály.
- subTask_Bounds.h - Alsó becslések előállításra szolgáló algoritmusok gyüjteménye
<br>

## Random példa generálás
A Sztandard inputról olvas be, három inputot:
- A példa tipus, egyenlőre mindig "TSP"
- Number of nodes (n): a generálandó gráf csúcsainak száma
- Edge ratio (mratio): a random generált élek és a csúcsok számának hányadosa, ha nagyobb mind az összes lehetséges él száma akkor $${n}\choose{2}$$
- Én az echo -e "TSP\n20\n10" | ./Rnd parancs paramétereinek modósításával szoktam futtatni

A generálás a következő módon zajlik:
- Veszünk n random egész pontot a síkon ezeknek egy új gráfban megfeleltetünk n csúcsot. Elkezdünk hozzáadni élek egy él hossza mindig az euklédeszi távolság felső egész része.
- Először Létrehozunk egy H - kört a gráfban, ez garantálja, az összefüggőséget, és, hogy legalább egy TSP megoldás létezzen
- A maradék éleket random generáljuk (nem generálunk multi éleket, ekkor újra sorsolunk)
- A Lemon graphwriter eszközeivel elmentjük a példát egy "digraph_tsp.lgf" nevű fájlba.

### AI use
A Kódot mivel mellékhatás nélküli egyszerű részfeladatot lát el, egy általam adott prompt és pszeudokód alapján a Google Gemini 1.5 Flash modell generálta én csak enyhe revíziókat, javítások végeztem, beépítettem a már létező programba.
<br>

## A solver osztály
A smallTSPsolver.cc része egy STSPsolver nevű template osztály ami a következő módon működik
- Az első template argument lehet Logging vagy Silent, ez dönti el, hogy a program írja a lépéseket egy "debug.log" fájlba, bármilyen beállítás esetén legfeljebb 1500 sort ír.
- A második template argument azt határozza meg milyen módon álapítsa meg a alsó korlátokat, jelenleg csak az "SST" opció működőképes.
- A harmadik opció $\alpha$-közelítő módszer ami már akkor is levág, ha a feladhoz tartozó alsó becslés $\alpha$-szorosa nagyobb az aktuális felső becslésnél. Mivel nem minden C++ verzió támogat double template argumentet ez a funkció constexpr segítségével van implementálva, az STSP::approx(double a) függvényt kell használni a közelítő hányados megadására.
- Az algoritmus argumentumként egyenlőre egy Lemon::Listtdigraph formátumú gráfot, rajta egy int tipusú 0, ... , n-1 csúcsindexelést egy nodemap fomrájában, és egy Arcmap<double> költségfüggvényt vár, opcionális egy kezdő felső korlát az optimumra (ha ez kissebb mint az optimum nem ad vissza megoldást).
- Megtalál egy 0-ból induló Hamilton utat, ami az ilyenek közt a költségfüggvényre nézve minimális.

A fő algoritmus egy Backtracking logikát követő bejárás, ami a maradék csúcsokon folyamatosan alsó becslést végez. Ha egy részfeladat már ki lett számolva, SparseMap<int, double> Lbounds{0.0} adattagban bitmaszk segítségével tárolja az alsó becslést, ennek a 0 alapértértelmezett értékét használom arra, hogy egy adott részfeladatról eldöntsem, már tartozik e hozzá alsó becslés, ez ha az alsó becslés sok 0-át ad jelenthet felesleges ellenörzést.
### Metódusok
- Konstruktor
A fentiekben leírt módon kapja példányosítja az osztált, mellet egyéb belső adattagokat is alaphelyzetbe állít
- doulbe solve() megoldja a feladatot
- void printroute(), vector<int> OPTroute(), double OPTval() visszadják, vagy printelik a várt módon a megoldást, értékét. Ha a solve előtt hívjuk exeptiont dobnak. Az út mindig -1 el kezdődik.

### AI use
A beolvasás a main függvényben, illetve a variadic template szintaktika a logfile írásban. Ezen kívűl semmi.
<br>

## Az alsó korlát header
Ez a fájl structokat tartalmaz amik, egy adott részfeladatra adnak alsó becslést. Fő célja az STSP solver segítése de próbáltam olyan absztraktcióval implementálni őket, hogy más feladatokra is a jövöben alkalmazhatók legyenek.
### SST (Shortest Spanning Tree)
Az alg a kruskal algoritmust használja minimális feszítőerdő megtalálására. Ez triviálisan jó alsó becslés
### Brute Force
Nem végez ellenörzést és a solverben használt constexpr miatt az egész constructbound lépés kioptimalizálódik.
### AP
  Az összefüggöségi követelmény elhagyásával, keres olyan részgráfot amiben, kettő csúcs (út eleje, vége) kivételével minden csúcs ki és befoka pontosan 1. Mivel a Lemon általam letöltött verziójában valamiért nincs MaxWeightedBipartiteMatching vagy ennek teljes párosítás kereső verziója. Az Edmonds munkáján alapuló MaxWeightedPerfectMatching osztályt használom.
### AI use
none
