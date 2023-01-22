#include "graf.h"

/*---------------------------------------------------------------------- Fa -----------------------------------------------------------------------------*/
// Utility
void Fa::level_elek_visit(int ind, double& osszeg) const
{
	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (this->is_leaf(szomszed.ind))
		{
			osszeg += szomszed.suly;
		}
		else
		{
			this->level_elek_visit(szomszed.ind, osszeg);
		}
	}
}

bool Fa::is_leaf(int node_ind) const
{
	/*
		Ellenorzi, ha egy csomopont level vagy sem. Egy csomopont akkor level, ha nincs egyetlen gyereke sem == szomszedsagi listaja ures.
	*/

	if (node_ind < 0 || node_ind >= this->n) throw WrongIndex();
	return (this->szomszedsagi_lista[node_ind].size() == 0);
}

// User functions
double Fa::level_elek_osszeg() const
{
	/*
		A level csomopontokhoz tartozo elek osszeget teriti.
	*/

	double osszeg = 0;

	this->level_elek_visit(this->gyoker, osszeg);

	return osszeg;
}

void Fa::torol_csomopont(int ind)
{
	// Ahhoz, hogy toroljunk egy csomopontot ki kell toroljunk minden el-t amivel csatlakozik
	auto it = this->el_lista.begin();
	while (it != this->el_lista.end())
	{
		if (it->kezd == ind || it->veg == ind)					// Toroljuk, ha csatlakozik
		{
			this->el_lista.erase(it);
		}
		else
		{
			it++;
		}
	}

	// Ha kell frissitsuk a csomopontok szamat akkor megtesszuk
	if (ind == this->n - 1)										// Ha utolso csomopont
	{
		this->n--;
	}

	this->el_lista_to_szomszedsagi_lista();
}

/*---------------------------------------------------------------------- Fa -----------------------------------------------------------------------------*/
/*---------------------------------------------------------------- Iranyitott Fa ------------------------------------------------------------------------*/
// Class
IranyitottFa::IranyitottFa()
{
	this->gyoker = 0;
}

IranyitottFa::IranyitottFa(string file_in)
{
	ifstream fin(file_in);

	if (fin.fail()) throw FileNotAccesible();

	fin >> this->n >> this->m >> this->gyoker;
	this->gyoker--;

	if ((this->n - 1) != this->m) throw GraphIsNotATree();

	El el;
	for (int i = 0; i < this->m; ++i)
	{
		fin >> el.kezd >> el.veg >> el.suly;
		el.kezd--; el.veg--;											// 0-tol valo indexeles
		this->el_lista.push_back(el);
	}

	this->el_lista_to_szomszedsagi_lista();

	if (this->van_kor()) throw GraphHasCycle();
}

IranyitottFa::IranyitottFa(const IranyitottFa& fa)
{
	this->n = fa.n; this->m = fa.m;
	this->el_lista = fa.el_lista;
	this->szomszedsagi_lista = fa.szomszedsagi_lista;
}

// Utility
void IranyitottFa::el_lista_to_szomszedsagi_lista()
{
	this->szomszedsagi_lista.resize(this->n);							// n darab csomopont
	for (list<SzListaElem>& lista : this->szomszedsagi_lista) lista.clear();

	for (const El& el : this->el_lista)
	{
		this->szomszedsagi_lista[el.kezd].push_back(SzListaElem(el.veg, el.suly));
	}
}

void IranyitottFa::LCA_DFS(int ind, int melyseg, vector<int>& node_depth, vector<int>& eulerian_tour, vector<int>& last, vector<bool>& latogatott) const
{
	latogatott[ind] = true;												// Standard DFS eljaras. Megjeloljuk, hogy itt mar jartunk
	eulerian_tour.push_back(ind);										// Az euler korbe helyezzuk a jelenlegi csomopontot
	node_depth.push_back(melyseg);										// A jelenlegi node melysege
	last[ind] = eulerian_tour.size() - 1;								// Megjegyezzuk az utolso indexet amin ez a node szerepel (Ha tobbszor szerepel akkor meg frissitjuk)

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			this->LCA_DFS(szomszed.ind, melyseg + 1, node_depth, eulerian_tour, last, latogatott);

			eulerian_tour.push_back(ind);									// Az euler korbe helyezzuk a jelenlegi csomopontot
			node_depth.push_back(melyseg);									// A jelenlegi node melysege
			last[ind] = eulerian_tour.size() - 1;							// Megjegyezzuk az utolso indexet amin ez a node szerepel (Ha tobbszor szerepel akkor meg frissitjuk)
		}
	}
}

// Operators
IranyitottFa& IranyitottFa::operator =(const IranyitottFa& fa)
{
	if (this != &fa)
	{
		this->el_lista = fa.el_lista;
		this->szomszedsagi_lista = fa.szomszedsagi_lista;
	}
	return *this;
}

ostream& operator <<(ostream& stream, const IranyitottFa& fa)
{
	return fa.kiir(stream);
}

// User functions
int IranyitottFa::center_node() const
{
	/*
		Egy csomopont akkor kozepso csomopont, ha kozepso vagy egyik kozepso csomopont a fa minden leghoszabb utjan.
		Egy fanak lehet egy vagy ketto kozepso csomopontja.
		Iranyitott grafot a jelenlegi implementacio iranyitatlankent kezeli.
	*/

	vector<int> csomopont_fokszam(this->n, 0), levelek;
	for (const El& el : this->el_lista)
	{
		csomopont_fokszam[el.kezd]++; csomopont_fokszam[el.veg]++;
	}

	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (csomopont_fokszam[csomopont] == 0 || csomopont_fokszam[csomopont] == 1)
		{
			levelek.push_back(csomopont);
		}
	}

	int level_darab = levelek.size();
	while (level_darab < this->n)
	{
		vector<int> uj_levelek = vector<int>();													// Egy ures vector a leveleknek
		for (const int& level_node : levelek)
		{
			for (const El& el : this->el_lista)
			{
				if (el.kezd == level_node)
				{
					if (csomopont_fokszam[el.veg]) csomopont_fokszam[el.veg]--;
					if (csomopont_fokszam[el.veg] == 1) uj_levelek.push_back(el.veg);			// Ha level lett
				}

				if (el.veg == level_node)
				{
					if (csomopont_fokszam[el.kezd]) csomopont_fokszam[el.kezd]--;
					if (csomopont_fokszam[el.kezd] == 1) uj_levelek.push_back(el.kezd);			// Ha level lett
				}
			}
			csomopont_fokszam[level_node] = 0;													// 'Toroljuk' a node-t
		}
		level_darab += uj_levelek.size();
		levelek = uj_levelek;
	}
	return levelek[0] + 1;
}

int IranyitottFa::legkisebb_kozos_os(int node_1, int node_2) const
{
	/*
		Legkisebb kozos os (LCA = lowest common ancestor) egy olyan csomopont amely ose mindket parameterkent
		adott csomopontnak es a tavolsaga a gyokertol a legnagyobb. (Tavolsaga a legnyaobb a gyokertol = legalacsonyabb os)

		Az algoritmus eloszor egy elofeldolgozast vegez amelyben letrehoz egy Euler kort es amely kor menten megjegyzi,
		hogy az adott csomopontoknak mekkora a melyseguk. Ezenkivul megjegyzi minden csomopontra, hogy mi az utolso
		elofordulasa. (Utolso elofordulas alatt azt ertjuk, hogy az Euler utvonalon mikor jelenik meg utoljara (last vector) )
		Ezutan pedig az adott csomopontok utolso indexei kozott a minimum index erteket teriti az algoritmus, amely az LCA.
	*/

	node_1--; node_2--;																// 0 -tol indexelunk
	if (node_1 < 0 || node_2 < 0 || node_1 >= this->n || node_2 >= this->n) throw WrongIndex();

	vector<bool> latogatott(this->n, false);
	vector<int> eulerian_tour, node_depth, last(this->n);

	this->LCA_DFS(this->gyoker, 0, node_depth, eulerian_tour, last, latogatott);	// Euler kor es csomopont melysegek kiszamitasa
	cout << node_depth;
	cout << eulerian_tour;
	cout << last;

	int min_depth = INT32_MAX, min_depth_ind = -1;									// Az algoritmus a legkisebb indexet teriti a ket csomopont altal meghatarozott intervallumbol
	int bal_ind = min(last[node_1], last[node_2]) , jobb_ind = max(last[node_1], last[node_2]);

	for (int ind = bal_ind; ind <= jobb_ind; ++ind)
	{
		if (node_depth[ind] < min_depth)
		{
			min_depth = node_depth[ind];
			min_depth_ind = eulerian_tour[ind];
		}
	}

	return min_depth_ind + 1;
}

/*---------------------------------------------------------------- Iranyitott Fa ------------------------------------------------------------------------*/
/*--------------------------------------------------------------- Iranyitatlan Fa -----------------------------------------------------------------------*/
// Class
IranyitatlanFa::IranyitatlanFa()
{
	this->gyoker = 0;
}

IranyitatlanFa::IranyitatlanFa(string file_in)
{
	ifstream fin(file_in);

	if (fin.fail()) throw FileNotAccesible();

	fin >> this->n >> this->m >> this->gyoker;
	this->gyoker--;

	if ((this->n - 1) != this->m) throw GraphIsNotATree();

	El el;
	for (int i = 0; i < this->m; ++i)
	{
		fin >> el.kezd >> el.veg >> el.suly;
		el.kezd--; el.veg--;											// 0-tol valo indexeles
		this->el_lista.push_back(el);
	}

	this->el_lista_to_szomszedsagi_lista();

	if (this->van_kor()) throw GraphHasCycle();
}

IranyitatlanFa::IranyitatlanFa(const IranyitatlanFa& fa)
{
	this->n = fa.n; this->m = fa.m;
	this->el_lista = fa.el_lista;
	this->szomszedsagi_lista = fa.szomszedsagi_lista;
}

// Operators
IranyitatlanFa& IranyitatlanFa::operator =(const IranyitatlanFa& fa)
{
	if (this != &fa)
	{
		this->el_lista = fa.el_lista;
		this->szomszedsagi_lista = fa.szomszedsagi_lista;
	}
	return *this;
}

ostream& operator <<(ostream& stream, const IranyitatlanFa& fa)
{
	return fa.kiir(stream);
}

// Utility
void IranyitatlanFa::el_lista_to_szomszedsagi_lista()
{
	this->szomszedsagi_lista.resize(this->n);							// n darab csomopont
	for (list<SzListaElem>& lista : this->szomszedsagi_lista) lista.clear();

	for (const El& el : this->el_lista)
	{
		this->szomszedsagi_lista[el.kezd].push_back(SzListaElem(el.veg, el.suly));
		this->szomszedsagi_lista[el.veg].push_back(SzListaElem(el.kezd, el.suly));
	}
}

// User functions
int IranyitatlanFa::center_node() const
{
	/*
		Egy csomopont akkor kozepso csomopont, ha kozepso vagy egyik kozepso csomopont a fa minden leghoszabb utjan.
		Egy fanak lehet egy vagy ketto kozepso csomopontja.
	*/

	vector<int> csomopont_fokszam(this->n, 0), levelek;
	for (const El& el : this->el_lista)
	{
		csomopont_fokszam[el.kezd]++; csomopont_fokszam[el.veg]++;
	}

	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (csomopont_fokszam[csomopont] == 0 || csomopont_fokszam[csomopont] == 1)
		{
			levelek.push_back(csomopont);
		}
	}

	int level_darab = levelek.size();
	while (level_darab < this->n)
	{
		vector<int> uj_levelek = vector<int>();													// Egy ures vector a leveleknek
		for (const int& level_node : levelek)
		{
			for (const SzListaElem& szomszed : this->szomszedsagi_lista[level_node])
			{
				if(csomopont_fokszam[szomszed.ind]) csomopont_fokszam[szomszed.ind]--;
				if (csomopont_fokszam[szomszed.ind] == 1) uj_levelek.push_back(szomszed.ind);	// Ha level lett
			}
			csomopont_fokszam[level_node] = 0;													// 'Toroljuk' a node-t
		}
		level_darab += uj_levelek.size();
		levelek = uj_levelek;
	}
	return levelek[0] + 1;
}
/*--------------------------------------------------------------- Iranyitatlan Fa -----------------------------------------------------------------------*/

/*------------------------------------------------------------------ Binaris Fa -------------------------------------------------------------------------*/
// Class
BinarisFa::BinarisFa()
	: IranyitottFa()
{

}

BinarisFa::BinarisFa(string file_in)
	: IranyitottFa(file_in)
{
	for (const list<SzListaElem>& szomszedok : this->szomszedsagi_lista)
	{
		if (szomszedok.size() > 2) throw TreeIsNotBinary();				// Akkor nem binaris fa, ha barmely csomopontnak tobb mint 2 gyereke van
	}
}

BinarisFa::BinarisFa(const vector<int>& Prufer_kod)
	: BinarisFa()
{
	this->n = Prufer_kod.size() + 1;									// A Prufer kodolas mindig n-1 hosszu, ahol n a csomopontok szama
	this->gyoker = Prufer_kod[Prufer_kod.size() - 1];					// A Prufer kodolasnal az utolso elem mindig a gyokeret jelkepezi

	vector<int> Prufer = Prufer_kod;									// Egy masolat amin tudunk modositani
	for (int& elem : Prufer) elem--;									// 0-tol indexelunk
	for (int iteracio = 0; iteracio < Prufer.size(); ++iteracio)
	{
		vector<bool> jelenleg_eleme(this->n, false);
		int legkisebb_nem_eleme = INT32_MAX;

		for (const int& elem : Prufer) jelenleg_eleme[elem] = true;

		// Megtalal legkisebb nem eleme index a jelenlegi sorozatbol
		for (int ind = 0; ind < jelenleg_eleme.size(); ++ind)
		{
			if (!jelenleg_eleme[ind])
			{
				legkisebb_nem_eleme = ind;
				break;
			}
		}

		// A vektor elso elemevel es megtalalt legkisebb index-el leterhozunk egy elet
		El uj_el(Prufer[0], legkisebb_nem_eleme, 1);
		*this = *this + uj_el;											// Hozzaadjuk az elet a jelenlegi objektumhoz

		// Most eltoljuk a vektort 1-el es a vegere helyezzuk a kapott legkisebb indexet
		for (int ind = 0; ind < Prufer.size() - 1; ++ind)
		{
			Prufer[ind] = Prufer[ind + 1];
		}
		Prufer[Prufer.size() - 1] = legkisebb_nem_eleme;
	}
}

BinarisFa::BinarisFa(const BinarisFa& fa)
	: IranyitottFa(fa)
{
}

// Operators
BinarisFa& BinarisFa::operator =(const BinarisFa& fa)
{
	if (this != &fa)
	{
		this->el_lista = fa.el_lista;
		this->szomszedsagi_lista = fa.szomszedsagi_lista;
	}
	return *this;
}

BinarisFa BinarisFa::operator +(const El& el)
{
	BinarisFa fa(*this);

	for (const El& jelen_el : fa.el_lista)
	{
		if (jelen_el == el) throw EdgeAlreadyExists();
	}

	fa.el_lista.push_back(el);
	fa.el_lista_to_szomszedsagi_lista();

	return fa;
}

BinarisFa& BinarisFa::operator +=(const El& el)
{
	return *this = *this + el;
}

ostream& operator <<(ostream& stream, const BinarisFa fa)
{
	return fa.kiir(stream);
}

// Utility
int BinarisFa::csomopont_magassag_rekurziv(int ind) const
{
	if (ind == -1) return 0;											// Ha level
	if (this->szomszedsagi_lista[ind].size() == 0) return 0;			// Ha level

	int szomszed_1 = -1, szomszed_2 = -1;
	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (szomszed_1 == -1) szomszed_1 = szomszed.ind;
		else szomszed_2 = szomszed.ind;
	}
	return max(this->csomopont_magassag_rekurziv(szomszed_1), this->csomopont_magassag_rekurziv(szomszed_2)) + 1;
}

void BinarisFa::preorder_visit(int ind, vector<int>& csomopontok) const
{
	csomopontok.push_back(ind);

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		this->preorder_visit(szomszed.ind, csomopontok);
	}
}

void BinarisFa::inorder_visit(int ind, vector<int>& csomopontok) const
{
	auto it = this->szomszedsagi_lista[ind].begin();

	if (this->szomszedsagi_lista[ind].size() > 0)
		this->inorder_visit(it->ind, csomopontok);

	csomopontok.push_back(ind);

	if (this->szomszedsagi_lista[ind].size() > 1)
	{
		++it;
		this->inorder_visit(it->ind, csomopontok);
	}
}

void BinarisFa::postorder_visit(int ind, vector<int>& csomopontok) const
{
	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		this->preorder_visit(szomszed.ind, csomopontok);
	}

	csomopontok.push_back(ind);
}

void BinarisFa::torol_level(int level_ind)
{
	if (!this->is_leaf(level_ind)) throw WrongIndex();

	auto it = this->el_lista.begin();
	while (it != this->el_lista.end())
	{
		if (it->kezd == level_ind || it->veg == level_ind)
		{
			this->el_lista.erase(it);
			it = this->el_lista.begin();
		}
		else it++;
	}

	this->m = el_lista.size();
	this->el_lista_to_szomszedsagi_lista();
}

// User functions
int BinarisFa::csomopont_magassag(int ind) const
{
	/*
		Egy csomopont magassagat teriti.
		Csomopont magassaga megegyezik a ket reszfajanak magassaganak maximumanak + 1. Rekurzivan definialva.
	*/

	ind--;																// 0-tol indexelunk
	if (ind < 0 || ind >= this->n) throw WrongIndex();
	return this->csomopont_magassag_rekurziv(ind);
}

int BinarisFa::get_parent(int ind) const
{
	for (const El& el : this->el_lista)
	{
		if (el.veg == ind) return el.kezd;
	}
	return -1;
}

vector<int> BinarisFa::preorder_bejaras() const
{
	vector<int> csomopontok;

	this->preorder_visit(this->gyoker, csomopontok);

	return csomopontok;
}

vector<int> BinarisFa::inorder_bejaras() const
{
	vector<int> csomopontok;

	this->inorder_visit(this->gyoker, csomopontok);

	return csomopontok;
}

vector<int> BinarisFa::postorder_bejaras() const
{
	vector<int> csomopontok;

	this->postorder_visit(this->gyoker, csomopontok);

	return csomopontok;
}

vector<int> BinarisFa::Prufer_kodolas() const
{
	/*
		A Prufer kodolas egy olyan kodolas, amely egy binaris fa teljes allapotat le tudja irni n - 1 szammal. (Ahol n a csomopontok szama)
		Az utolso szam a kapott kodolasban mindig a fa gyokere.
	*/

	BinarisFa fa(*this);
	int csomopont_db = fa.n;
	vector<int> kodolas;
	vector<bool> hasznalt_csomopont(csomopont_db, false);

	for (int level = 0; level < csomopont_db - 1; ++level)
	{
		int legkisebb_cimkeju_level = INT32_MAX;
		for (int csomopont = 0; csomopont < csomopont_db; ++csomopont)
		{
			if (fa.is_leaf(csomopont) && !hasznalt_csomopont[csomopont] && csomopont < legkisebb_cimkeju_level)
			{
				legkisebb_cimkeju_level = csomopont;
			}
		}

		// Behelyezzuk a level szulojet
		kodolas.push_back(fa.get_parent(legkisebb_cimkeju_level));
		hasznalt_csomopont[legkisebb_cimkeju_level] = true;
		
		// Toroljuk a levelet
		fa.torol_level(legkisebb_cimkeju_level);
	}

	return kodolas;
}
/*------------------------------------------------------------------ Binaris Fa -------------------------------------------------------------------------*/
