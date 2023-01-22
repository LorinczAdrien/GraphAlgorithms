#pragma once

#include <iostream>
#include <fstream>

#include <utility>
#include <tuple>

#include <vector>
#include <list>
#include <string>
#include <queue>
#include <stack>

using namespace std;

// Osztaly deklaraciok
class GridGraph;

struct El
{
	int kezd, veg;
	double suly;

	El(int kezd = 0, int veg = 0, double suly = 0.0) { this->kezd = kezd; this->veg = veg; this->suly = suly; }
	bool operator ==(const El& el) const { return ((this->kezd == el.kezd) && (this->veg == el.veg) && (this->suly == el.suly)); }
	bool operator !=(const El& el) const { return ((this->kezd != el.kezd) || (this->veg != el.veg) || (this->suly != el.suly)); }
};

struct SzListaElem
{
	int ind;
	double suly;

	SzListaElem(int szomszed = -1, double suly = 0.0) { this->ind = szomszed; this->suly = suly; }
	bool operator ==(const SzListaElem& szomszed) const { return this->ind == szomszed.ind && this->suly == szomszed.suly; }
	bool operator !=(const SzListaElem& szomszed) const { return this->ind != szomszed.ind || this->suly != szomszed.suly; }
};

struct PrimNode
{
	double distance;
	int parent;
	bool ismeretlen;

	PrimNode(double distance = DBL_MAX, int parent = -1, bool ismeretlen = true) { this->distance = distance; this->parent = parent; this->ismeretlen = ismeretlen; }
};

struct DijkstraNode
{
	int csomopont;
	double tavolsag;

	DijkstraNode(int csomopont = -1, double tavolsag = 0.0) { this->csomopont = csomopont, this->tavolsag = tavolsag; }
};

struct PumpaloNode
{
	int magassag;
	double tobbletfolyam;

	PumpaloNode(int magassag = 0, double ertek = 0) { this->magassag = magassag; this->tobbletfolyam = ertek; }
};

struct TevekenysegNode
{
	double vegrehajtasi_ido;
	double legkorabbi_kezdes, legkesobbi_kezdes;									// Legkorabbi es legkesobbi idopont amikor a tevekenyseg megkezdodhet
	double legkorabbi_befejezodes, legkesobbi_befejezodes;							// Legkorabbi es legkesobbi idopont amikor a tevekenyseg befejezodhet

	TevekenysegNode(double ido = 0, double korabbi_kezdes = 0, double kesobbi_kezdes = 0, double korabbi_befejezodes = DBL_MAX, double kesobbi_befejezodes = DBL_MAX)
	{
		this->vegrehajtasi_ido = ido; this->legkorabbi_kezdes = korabbi_kezdes; this->legkesobbi_kezdes = kesobbi_kezdes; this->legkorabbi_befejezodes = korabbi_befejezodes; this->legkesobbi_befejezodes = kesobbi_befejezodes;
	}
};

struct Subset
{
	int parent, rank;
	Subset(int parent, int rank) { this->parent = parent; this->rank = rank; }
};

class Graf
{
protected:
	int n, m;
	list<El> el_lista;
	vector< list<SzListaElem> > szomszedsagi_lista;
	vector< vector<double> > szomszedsagi_matrix, incidencia_matrix;

	virtual void el_lista_to_szomszedsagi_matrix() {}
	virtual void szomszedsagi_matrix_to_incidencia_matrix() {}
	virtual void incidencia_matrix_to_szomszedsagi_lista() {}
	virtual void szomszedsagi_lista_to_el_lista() {}

	virtual void torol_csomopont(int ind);

	virtual bool melysegi_visit_kor(int ind, vector<bool>& latogatott, vector<bool>& jelenlegi_ut) const;
	void topologiai_visit(int ind, vector<bool>& latogatott, vector<int>& topo) const;
	void melysegi_visit(int ind, vector<bool>& latogatott, vector<int>& bejart) const;
	void nodes_to_node_recursive(const vector<int>& parent, vector<int>& ut, int ind) const;
	vector<int> find_path_of_nodes_to_node(const vector<int>& parent, int ind) const;
	void megkeres_utvonal(int ind, int veg, list<El>& utvonal, const vector<int>& parent) const;

public:
	// Kivetelek
	class FileNotAccesible {};
	class WrongIndex {};
	class WrongValue {};
	class GraphHasCycle {};
	class GraphHasNegativeCycle {};
	class EdgeDoesNotExist {};
	class EdgeAlreadyExists {};

	// Graf tulajdonsagok
	bool regularis() const;																			// TRUE, ha a graf regularis, ellenkezo esetben FALSE
	virtual bool van_kor() const;																	// TRUE, ha a grafban van kor, ellenkezo esetben FALSE
	bool van_negativ_kor(vector<double> tavolsag = vector<double>()) const;							// TRUE, ha a grafban van negativ kor, ellenkezo esetben FALSE

	// Bejarasok
	vector<int> szelessegi_bejaras(int indulas_ind = 0) const;										// Szelessegi bejaras, amely visszateriti a bejart csomopontokat (Alapertelmezetten minden csomopontbol bejarja a grafot)
	vector<int> melysegi_bejaras(int indulas_ind = 0) const;										// Melysegi bejaras, amely visszateriti a bejart csomopontokat (Alapertelmezetten minden csomopontbol bejarja a grafot)
	
	// Specialis csomopontok, stb...
	vector<int> vegpontok() const;																	// Visszateriti a grafban talalhato vegpontokat
	vector<int> izolalt_csomopontok() const;														// Visszateriti a grafban talalhato izolalt csomopontokat
	vector<int> k_rendu_ismerosok(int ind, int k) const;											// Visszateriti egy adott csomopontbol k tavolsagra levo csomopontokat
	vector<int> topologiai_rendezes() const;														// Visszateriti a topologiailag rendezett csomopontokat

	// Legrovidebb ut es hossza egy adott csomopont es minden tobbi csomopont kozott
	pair< vector<double>, vector<int> > Moore_SP(int indulas_ind) const;							// Visszaterit egy minimalis tavolsag es szulo vector-t (A grafot sulyozatlan grafkent kezeli -> Minden el egy standard meretet hataroz meg)
	pair< vector<double>, vector<int> > Dijkstra_SP(int indulas_ind) const;							// Visszaterit egy minimalis tavolsag es szulo vector-t (Nem mukodik negativ el sulyok eseten)
	pair< vector<double>, vector<int> > Bellman_Ford_SP(int indulas_ind) const;						// Visszaterit egy minimalis tavolsag es szulo vector-t (Ha van negativ kor akkor kivetelt valt ki: GraphHasNegativeCycle)

	void kiir_el_lista() const;
	void kiir_szomszedsagi_matrix() const;
	void kiir_incidencia_matrix() const;
	void kiir_szomszedsagi_lista() const;
	ostream& kiir(ostream& stream) const;															// cout << -hoz 
};

class IranyitatlanGraf : public virtual Graf
{
protected:
	void el_lista_to_szomszedsagi_matrix();
	void szomszedsagi_matrix_to_incidencia_matrix();
	void incidencia_matrix_to_szomszedsagi_lista();
	void szomszedsagi_lista_to_el_lista();

	// DFS/BFS segedek
	bool melysegi_visit_kor(int volt_index, int ind, vector<bool>& latogatott, vector<bool>& jelenlegi_ut) const;
	void melysegi_visit_Tarjan_BiConnect(int jelen_ind, int szulo_ind, int& time, vector<int>& id, vector<int>& low, stack<El>& verem, vector< list<El> >& komponensek, list<int>& elvago_pontok, list<El>& hidak) const;
	int get_min_Prim(const vector<PrimNode>& csomopontok) const;									// Minimum csomopontba vezeto utat talalja meg (Prim algoritmus)

	void melysegi_visit_TSP_2_approximate(int ind, vector<bool>& latogatott, vector<int>& csomopontok) const;
	vector<int> melysegi_bejaras_TSP_2_approximate(int indulas_ind) const;

	// Union-Find segedek
	int find_union(int ind, vector<Subset>& subsets) const;
	void union_union(int x, int y, vector<Subset>& subsets) const;

public:
	IranyitatlanGraf();																				// Alapertelmezett konstruktor
	IranyitatlanGraf(string file_name);																// Beolvasassal inicializalo konstruktor
	IranyitatlanGraf(const IranyitatlanGraf& graf);													// Masolo konstruktor
	IranyitatlanGraf(const list<El>& elek);															// El listabol letrehozo konstruktor
	IranyitatlanGraf(const GridGraph& grid_graph);													// Grid grafbol iranyitalan grafot letrehozo konstruktor

	// Operator tulterhelesek
	IranyitatlanGraf& operator =(const IranyitatlanGraf& graf);										// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)
	
	IranyitatlanGraf operator +(const El& el);														// + operator tulterhelese, amely egy elet helyez a grafba, ha mar nem letezik (Exception: EdgeAlreadyExists) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitatlanGraf operator +(const list<El>& elek);												// + operator tulterhelese, amely a parameterkent megadott el listaban levo eleket behelyezi a grafba, ha meg nem leteznek
	IranyitatlanGraf operator +(const IranyitatlanGraf& graf);										// + operator tulterhelese, amely a parameterkent megadott grafban levo eleket behelyezi a grafba, ha meg nem leteznek
	IranyitatlanGraf& operator +=(const El& el);													// += operator tulterhelese, amely egy elet helyez a jelenlegi graba, ha mar nem letezik (Exception: EdgeAlreadyExists) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitatlanGraf& operator +=(const list<El>& elek);											// += operator tulterhelese, amely a parameterkent megadott el listaban levo eleket behelyezi a jelenlegi grafba, ha meg nem leteznek
	IranyitatlanGraf& operator +=(const IranyitatlanGraf& graf);									// += operator tulterhelese, amely a parameterkent megadott grafban levo eleket behelyezi a jelenlegi grafba, ha meg nem leteznek
	
	IranyitatlanGraf operator -(const El& el);														// - operator tulterhelese, amely egy adott elet torol a grafbol, ha letezik (Exception ha nem letezik: EdgeDoesNotExist) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitatlanGraf operator -(const list<El>& elek);												// - operator tulterhelese, amely a parameterkent megadott el listaban levo eleket torli a grafban, ha leteznek
	IranyitatlanGraf operator -(const IranyitatlanGraf& graf);										// - operator tulterhelese, amely a parameterkent megadott grafban levo eleket torli a grafban, ha leteznek
	IranyitatlanGraf& operator -=(const El& el);													// -= operator tulterhelese, amely egy adott elet torol a jelenlegi grafbol, ha letezik (Exception ha nem letezik: EdgeDoesNotExist) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitatlanGraf& operator -=(const list<El>& graf);											// -= operator tulterhelese, amely a parameterkent megadott el listaban levo eleket torli a jelenlegi grafbol, ha leteznek
	IranyitatlanGraf& operator -=(const IranyitatlanGraf& graf);									// -= operator tulterhelese, amely a parameterkent megadott grafban levo eleket torli a jelenlegi grafbol, ha leteznek
	
	bool van_kor() const override;
	bool van_kor_union_find() const;

	// Algoritmusok

	// Minimalis feszitofa es koltsege
	pair< list<El>, double > Prim_MST() const;														// Minimalis feszitofa eleit es koltseget teriti vissza
	pair< list<El>, double > Kruskal_MST() const;													// Minimalis feszitofa eleit es koltseget teriti vissza
	pair< list<El>, double > Boruvka_MST() const;													// Minimalis feszitofa eleit es koltseget teriti vissza
	pair< list<El>, double > Forditott_torles_MST() const;											// Minimalis feszitofa eleit es koltseget teriti vissza

	// Ketszeresen osszefuggo komponensek, elvago pontok es hidak
	tuple< vector< list<El> >, list<int>, list<El> > Tarjan_BiConnect() const;						// Osszefuggo komponenseket, elvago pontokat, illetve hidakat terit

	// Traveling salesman problem
	pair< list<El>, double > TSP_2_approximate() const;												// Traveling salesman problem-ra egy legroszabb esetben 2x megkozelitest ad. Eleket es Hamilton kor hosszat terit
	
	// Legrovidebb ut es hossza egy sulyozatlan grafban ket csomopont kozott
	pair< list<El>, int > BFS_SP(int kezd, int veg) const;											// Legrovidebb utat meghatarozza ket csomopont kozott, sulyozatlan grafok eseten mukodik csak. Az utat es hosszat teriti
};

class IranyitottGraf : public virtual Graf
{
protected:
	// Atalakitasok
	void el_lista_to_szomszedsagi_matrix();
	void el_lista_to_szomszedsagi_lista();
	void szomszedsagi_matrix_to_incidencia_matrix();
	void incidencia_matrix_to_szomszedsagi_lista();
	bool letezik_szomszed(int ind, const SzListaElem& szomszed) const;
	void szomszedsagi_lista_to_el_lista();

	// Seged algoritmusoks
	void transzponal_el_lista();
	int get_elek_szama_from_szomszedsagi_matrix() const;

	// Melysegi es szelessegi bejaras segedek
	void melysegi_visit_plusz_minusz(int ind, vector<bool>& plusz, vector<bool>& minusz, vector<bool>& latogatott, bool pluszt_keres) const;
	void melysegi_visit_Kosaraju_verem(int ind, vector<bool>& latogatott, stack<int>& verem) const;
	void melysegi_visit_Kosaraju_komponensek(int ind, vector<bool>& latogatott, vector< list<int> >& komponensek) const;
	void melysegi_visit_Tarjan(int jelen_ind, int& time, vector< list<int> >& komponensek, stack<int>& verem, vector<int>& id, vector<int>& low, vector<bool>& vermen) const;
	bool szelessegi_bejaras_Edmonds(int forras, int nyelo, vector<int>& parent, vector<vector<double>>& rezidualis_graf);
	vector<int> szelessegi_bejaras_Edmonds_csomopontok(int indulas_ind, vector< vector<double> >& rezidualis_graf);
	
	// Floyd Warshall
	void init_Floyd_Warshall(vector< vector<double> >& tavolsag, vector< vector<int> >& parent) const;

	// Pumpalo algoritmus segedek
	bool emelheto_csomopont(int csomopont, int forras, int nyelo, const vector<PumpaloNode>& csomopontok, const vector< vector<double> >& rezidualis_graf);
	bool pumpalhato_el(const El& el, int forras, int nyelo, const vector<PumpaloNode>& csomopontok, const vector< vector<double> >& rezidualis_graf);
	void Pumpalo_init(int forras, vector<PumpaloNode>& csomopontok, vector< vector<double> >& rezidualis_graf);
	void Pumpalas(const El& el, vector<PumpaloNode>& csomopontok, vector< vector<double> >& rezidualis_graf);
	void Emeles(int csomopont, vector<PumpaloNode>& csomopontok, vector< vector<double> >& rezidualis_graf);

public:
	IranyitottGraf();																				// Alapertelmezett konstruktor
	IranyitottGraf(string file_name);																// Beolvasassal inicializalo konstruktor
	IranyitottGraf(const vector< vector<double> >& matrix);											// Szomszedsagi/Incidencia matrixxal inicializalo konstruktor
	IranyitottGraf(const IranyitottGraf& graf);														// Masolo konstruktor
	
	// Operator tulterhelesek
	IranyitottGraf& operator =(const IranyitottGraf& graf);											// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)

	IranyitottGraf operator +(const El& el);														// + operator tulterhelese, amely egy elet helyez a grafba, ha mar nem letezik (Exception: EdgeAlreadyExists) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitottGraf operator +(const list<El>& elek);												// + operator tulterhelese, amely a parameterkent megadott el listaban levo eleket behelyezi a grafba, ha meg nem leteznek
	IranyitottGraf operator +(const IranyitottGraf& graf);											// + operator tulterhelese, amely a parameterkent megadott grafban levo eleket behelyezi a grafba, ha meg nem leteznek
	IranyitottGraf& operator +=(const El& el);														// += operator tulterhelese, amely egy elet helyez a jelenlegi graba, ha mar nem letezik (Exception: EdgeAlreadyExists) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitottGraf& operator +=(const list<El>& elek);												// += operator tulterhelese, amely a parameterkent megadott el listaban levo eleket behelyezi a jelenlegi grafba, ha meg nem leteznek
	IranyitottGraf& operator +=(const IranyitottGraf& graf);										// += operator tulterhelese, amely a parameterkent megadott grafban levo eleket behelyezi a jelenlegi grafba, ha meg nem leteznek

	IranyitottGraf operator -(const El& el);														// - operator tulterhelese, amely egy adott elet torol a grafbol, ha letezik (Exception ha nem letezik: EdgeDoesNotExist) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitottGraf operator -(const list<El>& elek);												// - operator tulterhelese, amely a parameterkent megadott el listaban levo eleket torli a grafban, ha leteznek
	IranyitottGraf operator -(const IranyitottGraf& graf);											// - operator tulterhelese, amely a parameterkent megadott grafban levo eleket torli a grafban, ha leteznek
	IranyitottGraf& operator -=(const El& el);														// -= operator tulterhelese, amely egy adott elet torol a jelenlegi grafbol, ha letezik (Exception ha nem letezik: EdgeDoesNotExist) es ujra letrehozza a szomszedsagi es incidencia matrixot, meg a szomszedsagi listat
	IranyitottGraf& operator -=(const list<El>& graf);												// -= operator tulterhelese, amely a parameterkent megadott el listaban levo eleket torli a jelenlegi grafbol, ha leteznek
	IranyitottGraf& operator -=(const IranyitottGraf& graf);										// -= operator tulterhelese, amely a parameterkent megadott grafban levo eleket torli a jelenlegi grafbol, ha leteznek

	// Algoritmusok

	// Tranzitiv lezaras
	vector< vector<bool> > tranzitiv_lezaras_matrix() const;										// Egy boolean matrixot terit, amely tartalmazza, hogy a jelenlegi objektum melyik csomopontjaibol melyibe lehet eljutni (Tranzitiv lezaras)				
	IranyitottGraf tranzitiv_lezaras() const;														// Egy iranyitott graf tipusu objektumot terit, amely tartalmazza a jelenlegi objektum tranzitiv lezarasat.

	// Erossen osszefuggo komponensek
	vector< list<int> > plusz_minusz_algoritmus() const;											// Erossen osszefuggo komponenseket, es azok szamat teriti vissza (komponens_db = vector<vector<int>>.size())
	vector< list<int> > Kosaraju_algoritmusa() const;												// Erossen osszefuggo komponenseket, es azok szamat teriti vissza (komponens_db = vector<vector<int>>.size())
	vector< list<int> > Tarjan_StrongConnect() const;												// Erossen osszefuggo komponenseket, es azok szamat teriti vissza (komponens_db = vector<vector<int>>.size())
	
	// Erossen osszefuggo komponensek, elvago pontok es hidak

	// Tavolsag minden csomopontbol minden csomopontba (Parent matrixok is)
	pair< vector< vector<double> >, vector< vector<int> > > Johnson();								// Visszateriti a tavolsag es szulo matrixot (Minden csomopontonbol minden csomopontba)
	pair< vector< vector<double> >, vector< vector<int> > > Floyd_Warshall() const;					// Visszateriti a tavolsag es szulo matrixot (Minden csomopontonbol minden csomopontba)
	
	// Maximalis folyam es minimalis vagat ket halmaza
	tuple < double, vector<int>, vector<int> > Edmonds_Karp(int forras, int nyelo);					// Maximalis folyam erteket es a ket minimalis vagatot tartalmazo halmazot teriti
	double Pumpalo_algoritmus(int forras, int nyelo);												// Maximalis folyam erteket teriti
	
	// Kritikus ut masodik modell
	vector<TevekenysegNode> Kritikus_ut_masodik_modell(const vector<double>& vegrehajtasi_idok);	// Minden csomopontra vonatkozo kritkus ut informaciot tartalmazo vector-t terit
};

class Fa : public virtual Graf
{
protected:
	int gyoker;

	virtual void el_lista_to_szomszedsagi_lista() = 0;
	void torol_csomopont(int ind);																	// Tulterhelt fuggveny: Torol egy csomopontot a grafbol

	void level_elek_visit(int ind, double& osszeg) const;

	bool is_leaf(int node_ind) const;																// True, ha level

public:
	class GraphIsNotATree {};

	// Algoritmusok
	double level_elek_osszeg() const;																// A levelekhez vezeto elek sulyanak az osszeget teriti
};

class IranyitottFa : public virtual Fa, public virtual IranyitottGraf
{
protected:
	void el_lista_to_szomszedsagi_lista();
	void LCA_DFS(int ind, int melyseg, vector<int>& node_melyseg, vector<int>& eulerian_tour, vector<int>& last, vector<bool>& latogatott) const;

public:
	IranyitottFa();																					// Alapertelmezett konstruktor
	IranyitottFa(string file_name);																	// Beolvasassal inicializalo konstruktor
	IranyitottFa(const IranyitottFa& graf);															// Masolo konstruktor
	IranyitottFa& operator =(const IranyitottFa& fa);												// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)

	// Algoritmusok
	int center_node() const;																		// Egy fa eseten megadja a 'kozepso' csomopontot. Egy csomopont akkor 'kozepso', ha a kozepso vagy az egyik kozepso csomopont minden egyes leghoszabb utan a faban.
	int legkisebb_kozos_os(int node_1, int node_2) const;											// LCA vagyis legkisebb kozos os-t teriti. LCA egy olyan csomopont amely mindket adott csomopontnak a szuloje es a legalacsonyabb, tehat a legtavolabb van a gyokertol.
};

class IranyitatlanFa : public virtual Fa, public virtual IranyitatlanGraf
{
protected:
	void el_lista_to_szomszedsagi_lista();

public:
	class GraphIsNotATree {};

	IranyitatlanFa();																			// Alapertelmezett konstruktor
	IranyitatlanFa(string file_name);															// Beolvasassal inicializalo konstruktor
	IranyitatlanFa(const IranyitatlanFa& graf);													// Masolo konstruktor
	IranyitatlanFa& operator =(const IranyitatlanFa& fa);										// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)

	// Algoritmusok
	int center_node() const;																	// Egy fa eseten megadja a 'kozepso' csomopontot. Egy csomopont akkor 'kozepso', ha a kozepso vagy az egyik kozepso csomopont minden egyes leghoszabb utan a faban.
};

class BinarisFa : public virtual IranyitottFa
{
protected:
	int csomopont_magassag_rekurziv(int ind) const;

	void preorder_visit(int ind, vector<int>& csomopontok) const;
	void inorder_visit(int ind, vector<int>& csomopontok) const;
	void postorder_visit(int ind, vector<int>& csomopontok) const;

	void torol_level(int level_ind);

public:
	class TreeIsNotBinary {};

	BinarisFa();																					// Alapertelmezett konstruktor
	BinarisFa(string file_name);																	// Beolvasassal inicializalo konstruktor
	BinarisFa(const vector<int>& Prufer);															// Prufer koddal inicializalo konstruktor
	BinarisFa(const BinarisFa& graf);																// Masolo konstruktor
	
	// Operatorok
	BinarisFa& operator =(const BinarisFa& fa);														// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)

	BinarisFa operator +(const El& el);
	BinarisFa& operator +=(const El& el);

	int get_parent(int ind) const;

	// Csomopont magassaga
	int csomopont_magassag(int ind) const;															// Egy adott csomopont magassag (Exception: WrongIndex, ha helytelen indexet adunk meg)

	// Binaris fa bejarasok
	vector<int> preorder_bejaras() const;
	vector<int> inorder_bejaras() const;
	vector<int> postorder_bejaras() const;

	// Prufer kodolas
	vector<int> Prufer_kodolas() const;																// Egy n - 1 hosszusagu vektort terit, amely tartalmazza a jelenlegi objektum Prufer kodolasat
};

struct Coord2D
{
	int x, y;

	// Konstruktor
	Coord2D(int x = 0, int y = 0) { this->x = x; this->y = y; }

	double tavolsag_ket_koordinata_kozott(const Coord2D& cord)
	{
		return sqrt( (cord.x - this->x) * (cord.x - this->x) + (cord.y - this->y) * (cord.y - this->y) );
	}

	// Operators
	Coord2D& operator =(const Coord2D& kord) { if (this != &kord) { this->x = kord.x; this->y = kord.y; } return *this; }

	Coord2D operator +(const Coord2D& kord) const { return Coord2D(this->x + kord.x, this->y + kord.y); }
	Coord2D& operator +=(const Coord2D& kord) { return *this = *this + kord; }
	
	Coord2D operator -(const Coord2D& kord) const { return Coord2D(this->x - kord.x, this->y - kord.y); }
	Coord2D& operator -=(const Coord2D& kord) { return *this = *this - kord; }

	bool operator ==(const Coord2D& kord) const { return ((this->x == kord.x) && (this->y == kord.y)); }
	bool operator !=(const Coord2D& kord) const { return ((this->x != kord.x) || (this->y != kord.y)); }
};

struct AStarNode
{
	Coord2D poz, parent_poz;
	double local_distance, heuristic_distance, total_distance;
	bool ismeretlen;

	AStarNode(const Coord2D& jelen_poz = Coord2D(0, 0), const Coord2D& parent_poz = Coord2D(0, 0), const Coord2D& start = Coord2D(0, 0), const Coord2D& end = Coord2D(0, 0), bool ismeretlen = true) {
		this->poz = jelen_poz;														// A node jelenlegi pozicioja
		this->parent_poz = parent_poz;												// Megjegyezzuk, hogy honnan jottunk
		this->local_distance = poz.tavolsag_ket_koordinata_kozott(start);			// Tavolsag a kezdeti node-tol
		this->heuristic_distance = poz.tavolsag_ket_koordinata_kozott(end);			// Becsult tavolsag a vegso node-ig
		this->total_distance = this->local_distance + this->heuristic_distance;		// Teljes koltseg (F cost)
		this->ismeretlen = ismeretlen; 
	}

	AStarNode& operator =(const AStarNode& node)
	{
		if (this != &node) { this->poz = node.poz; this->parent_poz = node.parent_poz; this->local_distance = node.local_distance; this->heuristic_distance = node.heuristic_distance; this->total_distance = node.total_distance; this->ismeretlen = node.ismeretlen; } return *this;
	}
};

class GridGraph
{
protected:
	int sor, oszlop;
	vector< vector<int> > grid;

	bool is_valid_coordinate(const Coord2D& kord) const;
	list<Coord2D> get_szomszedok(const Coord2D& jelen_poz, int mode) const;							// Egy adott pozicioval szomszedos koordinatakat terit
	list<Coord2D> get_szabad_szomszedok(const Coord2D& jelen_poz, int mode) const;					// Egy adott pozicioval szomszedos koordinatakat terit, ahol az ertek 0
	void megkeres_ut(Coord2D jelen_poz, list<Coord2D>& ut, const Coord2D& keresett, const vector< vector<Coord2D> >& parent) const;

	Coord2D D1_to_D2(int poz) const;
	int D2_to_D1(const Coord2D& poz) const;

public:
	class FileNotAccesible {};
	class IndexOutOfBounds {};
	class WallAtCoordinate {};

	GridGraph(string file_in);																		// Beolvasassal inicializalo konstruktor
	GridGraph(int sor = 0, int oszlop = 0, const vector< vector<int> >& adott_grid = vector< vector<int> > (0));
	GridGraph& operator =(const GridGraph& grid_graf);												// Ertekado operator tulterhelese (A parameterkent adott graf teljes tartalmat belemasolja a jelenlegi objektumba)

	void set_ertek(int sor, int oszlop, int uj_ertek);
	int get_ertek(int sor, int oszlop) const;

	const vector< vector<int> >& get_grid() const { return this->grid; }							// Konstans referencia a graf gridjehez
	list<int> get_szabad_szomszedok_1D(int poz, int mode) const;									// Iranyitatlan grafba valo atvaltoztatashoz

	// Algoritmusok
	pair< list<Coord2D>, int > BFS_SP(const Coord2D& start, const Coord2D& end) const;				// A legrovidebb utat teriti ket koordinata kozott es ennek az utnak a hosszat
	pair< list<Coord2D>, int > A_star_SP(const Coord2D& start, const Coord2D& end) const;			// A legrovidebb utat teriti ket koordinata kozott es ennek az utnak a hosszat

	ostream& kiir(ostream& stream) const;															// << operator
};

// << operator
ostream& operator <<(ostream& stream, const IranyitatlanGraf& graf);
ostream& operator <<(ostream& stream, const IranyitottGraf& graf);
ostream& operator <<(ostream& stream, const GridGraph& graf);
ostream& operator <<(ostream& stream, const IranyitottFa& fa);
ostream& operator <<(ostream& stream, const IranyitatlanFa& fa);

// bonusz: << operator tulterhelesek
ostream& operator <<(ostream& stream, const list<El>& elek);
ostream& operator <<(ostream& stream, const list<Coord2D>& koordinatak);
ostream& operator <<(ostream& stream, const vector<int>& v);
ostream& operator <<(ostream& stream, const list<int>& v);
ostream& operator <<(ostream& stream, const vector<double>& v);
ostream& operator <<(ostream& stream, const vector< vector<bool> >& v);
ostream& operator <<(ostream& stream, const vector< vector<int> >& v);
ostream& operator <<(ostream& stream, const vector< list<int> >& v);
ostream& operator <<(ostream& stream, const vector< vector<double> >& v);
ostream& operator <<(ostream& stream, const vector< vector<Coord2D> >& v);