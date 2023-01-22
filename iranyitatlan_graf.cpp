#include "graf.h"

/*--------------------------------------------------------------- Iranyitatlan Graf -------------------------------------------------------------------------*/
// Class
IranyitatlanGraf::IranyitatlanGraf()
{
	this->n = this->m = 0;
}

IranyitatlanGraf::IranyitatlanGraf(string file_name)
{
	ifstream fin(file_name);

	if (fin.fail()) throw FileNotAccesible();							// Ha nem sikerult megnyitni a fajl-t

	fin >> this->n >> this->m;

	El el;
	for (int i = 0; i < this->m; ++i)
	{
		fin >> el.kezd >> el.veg >> el.suly;
		el.kezd--; el.veg--;											// 0-tol valo indexeles
		this->el_lista.push_back(el);
	}

	this->el_lista_to_szomszedsagi_matrix();
	this->szomszedsagi_matrix_to_incidencia_matrix();
	this->incidencia_matrix_to_szomszedsagi_lista();
}

IranyitatlanGraf::IranyitatlanGraf(const IranyitatlanGraf& graf)
{
	this->n = graf.n; this->m = graf.m;

	// El lista masolas
	this->el_lista = graf.el_lista;

	// Szomszedsagi lista masolas
	this->szomszedsagi_matrix.resize(this->n);
	for (int i = 0; i < this->n; ++i)
	{
		this->szomszedsagi_matrix[i].resize(this->n);
		for (int j = 0; j < this->n; ++j) this->szomszedsagi_matrix[i][j] = graf.szomszedsagi_matrix[i][j];
	}

	// Incidencia matrix masolas
	this->incidencia_matrix.resize(this->n);
	for (int i = 0; i < this->n; ++i)
	{
		this->incidencia_matrix[i].resize(this->m);
		for (int j = 0; j < this->m; ++j) this->incidencia_matrix[i][j] = graf.incidencia_matrix[i][j];
	}

	// Szomszedsagi list masolas
	this->szomszedsagi_lista.resize(this->n);
	for (int i = 0; i < this->n; ++i)
	{
		this->szomszedsagi_lista[i] = graf.szomszedsagi_lista[i];
	}

}

IranyitatlanGraf::IranyitatlanGraf(const list<El>& elek)
{
	this->el_lista = elek;
	this->m = static_cast<int>(this->el_lista.size());

	this->n = 0;
	// Megnezzuk, hogy n minimalisan mennyi kell legyen (a legnagyobb kezd/veg pont az el_lista-bol)
	for(const El& el : this->el_lista)
	{
		if (el.kezd + 1 > this->n) this->n = el.kezd + 1;
		if (el.veg + 1 > this->n) this->n = el.veg + 1;
	}

	this->el_lista_to_szomszedsagi_matrix();
	this->szomszedsagi_matrix_to_incidencia_matrix();
	this->incidencia_matrix_to_szomszedsagi_lista();
}

IranyitatlanGraf::IranyitatlanGraf(const GridGraph& grid_graph)
{
	const vector< vector<int> >& grid = grid_graph.get_grid();

	int csomopont_db = 0;
	for (const vector<int>& csomopontok : grid)					// Meghatarozzuk a csomopontok szamat
	{
		csomopont_db += csomopontok.size();
	}
	this->szomszedsagi_lista.resize(csomopont_db);

	for (int i = 0; i < grid.size(); ++i)
	{
		for (int j = 0; j < grid[i].size(); ++j)
		{
			int jelen_csomopont = this->n;
			list<int> szomszedok = grid_graph.get_szabad_szomszedok_1D(jelen_csomopont, 4);		// 4 - ha diagonal is akkor 8
			this->n++;											// Minden cella egy csomopont

			for (const int& szomszed : szomszedok)				// Minden egyes szomszed eseten behelyezzuk az eleket a szomszedsagi listaba
			{
				this->szomszedsagi_lista[jelen_csomopont].push_back(SzListaElem(szomszed, 1));
			}
		}
	}

	this->szomszedsagi_lista_to_el_lista();
	this->el_lista_to_szomszedsagi_matrix();
	this->szomszedsagi_matrix_to_incidencia_matrix();
}

// Operators =, +, +=, -, -=
IranyitatlanGraf& IranyitatlanGraf::operator =(const IranyitatlanGraf& graf)
{
	if (this != &graf)
	{
		this->n = graf.n; this->m = graf.m;

		this->el_lista = graf.el_lista;

		this->szomszedsagi_matrix.resize(n);
		for (int i = 0; i < this->n; ++i)
		{
			this->szomszedsagi_matrix[i].resize(n);
			for (int j = 0; j < this->n; ++j) this->szomszedsagi_matrix[i][j] = graf.szomszedsagi_matrix[i][j];
		}

		this->incidencia_matrix.resize(n);
		for (int i = 0; i < this->n; ++i)
		{
			this->incidencia_matrix[i].resize(m);
			for (int j = 0; j < this->m; ++j) this->incidencia_matrix[i][j] = graf.incidencia_matrix[i][j];
		}

		this->szomszedsagi_lista.resize(n);
		for (int i = 0; i < this->n; ++i)
		{
			this->szomszedsagi_lista[i] = graf.szomszedsagi_lista[i];
		}
	}

	return *this;
}

IranyitatlanGraf IranyitatlanGraf::operator +(const El& el)
{
	IranyitatlanGraf osszeadott(*this);

	for (const El& jelen_el : osszeadott.el_lista)
	{
		if (jelen_el == el) throw EdgeAlreadyExists();					// Ha mar letezik a megadott el
	}

	if (el.kezd + 1 > osszeadott.n) osszeadott.n = el.kezd + 1;
	if (el.veg + 1 > osszeadott.n) osszeadott.n = el.veg + 1;

	osszeadott.el_lista.push_back(el);

	osszeadott.m = static_cast<int>(osszeadott.el_lista.size());
	osszeadott.el_lista_to_szomszedsagi_matrix();
	osszeadott.szomszedsagi_matrix_to_incidencia_matrix();
	osszeadott.incidencia_matrix_to_szomszedsagi_lista();

	return osszeadott;
}

IranyitatlanGraf IranyitatlanGraf::operator +(const list<El>& elek)
{
	IranyitatlanGraf osszeadott(*this);

	for (const El& el : elek)											// Ha az adott grafban tobb csomopont van, akkor noveljuk a szomszedsagi matrix meretet es a csomopontok szamat
	{
		if (el.kezd + 1 > osszeadott.n) osszeadott.n = el.kezd + 1;
		if (el.veg + 1 > osszeadott.n) osszeadott.n = el.veg + 1;
	}
	osszeadott.szomszedsagi_matrix.resize(osszeadott.n);
	for (int i = 0; i < osszeadott.n; ++i) osszeadott.szomszedsagi_matrix[i].resize(osszeadott.n);

	for (const El& el : elek)
	{
		if (osszeadott.szomszedsagi_matrix[el.kezd][el.veg] == 0)			// Ha a jelenlegi objektumban nem letezik ott el, viszont a parameterkent adott objektumban igen, akkor 'behelyezzuk' ide is
		{
			osszeadott.el_lista.push_back(el);
			osszeadott.m = static_cast<int>(osszeadott.el_lista.size());

			osszeadott.szomszedsagi_matrix[el.kezd][el.veg] = el.suly;
			osszeadott.szomszedsagi_matrix[el.veg][el.kezd] = el.suly;
		}
	}

	osszeadott.szomszedsagi_matrix_to_incidencia_matrix();
	osszeadott.incidencia_matrix_to_szomszedsagi_lista();

	return osszeadott;
}

IranyitatlanGraf IranyitatlanGraf::operator +(const IranyitatlanGraf& graf)
{
	IranyitatlanGraf osszeadott(*this);
	return osszeadott + graf.el_lista;
}

IranyitatlanGraf& IranyitatlanGraf::operator +=(const El& el)
{
	return *this = *this + el;
}

IranyitatlanGraf& IranyitatlanGraf::operator +=(const list<El>& elek)
{
	return *this = *this + elek;
}

IranyitatlanGraf& IranyitatlanGraf::operator +=(const IranyitatlanGraf& graf)
{
	return *this = *this + graf;
}

IranyitatlanGraf IranyitatlanGraf::operator -(const El& el)
{
	IranyitatlanGraf kivont(*this);

	bool edge_exists = false;
	for (list<El>::const_iterator it = kivont.el_lista.begin(); it != kivont.el_lista.end(); ++it)
	{
		if (*it == el)										// El == operator tul van terhelve
		{
			kivont.el_lista.erase(it);						// Kitoroljuk az elet
			edge_exists = true;
			break;
		}
	}

	if (!edge_exists) throw EdgeDoesNotExist();

	kivont.m = static_cast<int>(kivont.el_lista.size());
	kivont.el_lista_to_szomszedsagi_matrix();

	kivont.szomszedsagi_matrix_to_incidencia_matrix();
	kivont.incidencia_matrix_to_szomszedsagi_lista();

	return kivont;
}

IranyitatlanGraf IranyitatlanGraf::operator -(const list<El>& elek)
{
	IranyitatlanGraf kivont(*this);

	for (const El& el : elek)
	{
		if ((el.kezd < kivont.n) && (el.veg < kivont.n) && kivont.szomszedsagi_matrix[el.kezd][el.veg])
		{
			kivont.szomszedsagi_matrix[el.kezd][el.veg] = 0;
			kivont.szomszedsagi_matrix[el.veg][el.kezd] = 0;

			int uj_n = 0;
			for (const El& kivont_el : kivont.el_lista)				// Megnezzuk ha csokkent a csomopontok szama (Persze lehetnek izolalt csomopontok is viszont igy ezt elkeruljuk)
			{
				if (kivont_el.kezd + 1 > uj_n) uj_n = kivont_el.kezd + 1;
				if (kivont_el.veg + 1 > uj_n) uj_n = kivont_el.veg + 1;
			}
			kivont.n = uj_n;
			kivont.m--;
		}
	}
	kivont.szomszedsagi_matrix.resize(kivont.n);

	kivont.szomszedsagi_matrix_to_incidencia_matrix();
	kivont.incidencia_matrix_to_szomszedsagi_lista();
	kivont.szomszedsagi_lista_to_el_lista();

	return kivont;
}

IranyitatlanGraf IranyitatlanGraf::operator -(const IranyitatlanGraf& graf)
{
	IranyitatlanGraf kivont(*this);
	return kivont - graf.el_lista;
}

IranyitatlanGraf& IranyitatlanGraf::operator -=(const El& el)
{
	return *this = *this - el;
}

IranyitatlanGraf& IranyitatlanGraf::operator -=(const list<El>& elek)
{
	return *this = *this - elek;
}

IranyitatlanGraf& IranyitatlanGraf::operator -=(const IranyitatlanGraf& graf)
{
	return *this = *this - graf;
}

// Utility
void IranyitatlanGraf::el_lista_to_szomszedsagi_matrix()
{
	this->szomszedsagi_matrix.resize(this->n);
	for (int i = 0; i < this->n; ++i)
	{
		this->szomszedsagi_matrix[i].resize(this->n);
		for (int j = 0; j < this->n; ++j)
		{
			this->szomszedsagi_matrix[i][j] = 0;
		}
	}

	for(const El& el : this->el_lista)
	{
		this->szomszedsagi_matrix[el.kezd][el.veg] = el.suly;
		this->szomszedsagi_matrix[el.veg][el.kezd] = el.suly;
	}
}

void IranyitatlanGraf::szomszedsagi_matrix_to_incidencia_matrix()
{
	if (this->incidencia_matrix.size() != this->n)
	{
		this->incidencia_matrix.resize(this->n);
		for (int i = 0; i < this->n; ++i)
		{
			this->incidencia_matrix[i].resize(this->m);
			for (int j = 0; j < incidencia_matrix[i].size(); ++j)
			{
				this->incidencia_matrix[i][j] = 0;
			}
		}
	}

	int jelen_el = 0;																// Szamoljuk, hogy hany elt talaltunk eddig
	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		for (int oszlop_ind = sor_ind + 1; oszlop_ind < this->n; ++oszlop_ind)		// sor_ind + 1, ugyanis eleg, ha a fo atlo felett nezzuk az ertekeket
		{
			if (this->szomszedsagi_matrix[sor_ind][oszlop_ind] != 0)
			{
				this->incidencia_matrix[sor_ind][jelen_el] = this->szomszedsagi_matrix[sor_ind][oszlop_ind];
				this->incidencia_matrix[oszlop_ind][jelen_el] = this->szomszedsagi_matrix[sor_ind][oszlop_ind];

				jelen_el++;															// Csak itt noveljuk, mivel 0-tol indexelunk
			}
		}
	}
}

void IranyitatlanGraf::incidencia_matrix_to_szomszedsagi_lista()
{
	if (this->szomszedsagi_lista.size() != this->n) this->szomszedsagi_lista.resize(this->n);					// n csomopont
	for (int i = 0; i < this->n; i++) this->szomszedsagi_lista[i].clear();

	for (int oszlop_ind = 0; oszlop_ind < this->m; ++oszlop_ind)
	{
		bool talalt = false;
		int sor_ind = 0, elso_ind = -1, masodik_ind = -1;

		while ((sor_ind < this->n) && !talalt)
		{
			if ((this->incidencia_matrix[sor_ind][oszlop_ind] != 0) && (elso_ind == -1))
			{
				elso_ind = sor_ind;
			}
			else
			{
				if ((this->incidencia_matrix[sor_ind][oszlop_ind] != 0) && (masodik_ind == -1))
				{
					masodik_ind = sor_ind;
					talalt = true;								// Tehat megtalaltuk, hogy ez az el melyik ket csomopontot kot ossze (megalhat a while)
				}
			}

			sor_ind++;
		}

		this->szomszedsagi_lista[elso_ind].push_back(SzListaElem(masodik_ind, this->incidencia_matrix[elso_ind][oszlop_ind]));	// Betesszuk a szomszedot es az el sulyat
		this->szomszedsagi_lista[masodik_ind].push_back(SzListaElem(elso_ind, this->incidencia_matrix[elso_ind][oszlop_ind]));
	}
}

void IranyitatlanGraf::szomszedsagi_lista_to_el_lista()
{
	this->el_lista = list<El>();

	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		for(const SzListaElem& szomszed : this->szomszedsagi_lista[sor_ind])
		{
			if (sor_ind < szomszed.ind)
			{
				// Behelyezzuk a csomopont indexet, szomszedjat es az el sulyat
				this->el_lista.push_back(El(sor_ind, szomszed.ind, szomszed.suly));
			}
		}
	}
	this->m = static_cast<int>(this->el_lista.size());
}

int IranyitatlanGraf::get_min_Prim(const vector<PrimNode>& csomopontok) const
{
	int minimum_ind = -1;
	double minimum_ertek = DBL_MAX;

	for (int csomopont_ind = 0; csomopont_ind < csomopontok.size(); ++csomopont_ind)
	{
		if (csomopontok[csomopont_ind].ismeretlen && csomopontok[csomopont_ind].distance < minimum_ertek)
		{
			minimum_ertek = csomopontok[csomopont_ind].distance;
			minimum_ind = csomopont_ind;
		}
	}

	return minimum_ind;
}

void IranyitatlanGraf::melysegi_visit_TSP_2_approximate(int ind, vector<bool>& latogatott, vector<int>& csomopontok) const
{
	latogatott[ind] = true;

	csomopontok.push_back(ind);
	bool level = true;
	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			level = false;
			this->melysegi_visit_TSP_2_approximate(szomszed.ind, latogatott, csomopontok);
		}
	}
	if (!level) csomopontok.push_back(ind);
}

vector<int> IranyitatlanGraf::melysegi_bejaras_TSP_2_approximate(int indulas_ind) const
{
	if (indulas_ind < 0 || indulas_ind > this->n) throw WrongIndex();

	vector<int> csomopontok;
	vector<bool> latogatott(this->n, false);

	indulas_ind--;
	this->melysegi_visit_TSP_2_approximate(indulas_ind, latogatott, csomopontok);

	return csomopontok;
}

bool IranyitatlanGraf::melysegi_visit_kor(int volt_ind, int ind, vector<bool>& latogatott, vector<bool>& jelenlegi_ut) const
{
	bool van_kor = false;
	latogatott[ind] = true; jelenlegi_ut[ind] = true;

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (szomszed.ind != volt_ind)								// Elkeruljuk az 1 hosszu 'korok' megtalalasat
		{
			if (jelenlegi_ut[szomszed.ind])							// Ha a jelenlegi utunkban mar szerepelt ez a csomopont (jelenleg egy szomszed), akkor talaltunk egy kort
			{
				return true;
			}
			else
			{
				if (!latogatott[szomszed.ind])
				{
					van_kor = this->melysegi_visit_kor(ind, szomszed.ind, latogatott, jelenlegi_ut);
					if (van_kor) return true;
				}
			}
		}
	}

	jelenlegi_ut[ind] = false;								// Visszalepeskor frissitjuk a jelenlegi utat (mint backtrack eseten)

	return van_kor;
}

void IranyitatlanGraf::melysegi_visit_Tarjan_BiConnect(int jelen_ind, int szulo_ind, int& time, vector<int>& id, vector<int>& low, stack<El>& verem, vector< list<El> >& komponensek, list<int>& elvago_pontok, list<El>& hidak) const
{
	id[jelen_ind] = low[jelen_ind] = ++time;

	int gyerekek_szama = 0;
	bool elvago_pont = false;

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_ind])
	{
		El jelen_el_to_szomszed(jelen_ind, szomszed.ind, this->szomszedsagi_matrix[jelen_ind][szomszed.ind]);

		if (id[szomszed.ind] == -1)																// Ha meg nem voltunk ott
		{
			gyerekek_szama++;																	// A gyerek szamat az elvago pontok meghatarozasanal hasznaljuk
			verem.push(jelen_el_to_szomszed);

			this->melysegi_visit_Tarjan_BiConnect(szomszed.ind, jelen_ind, time, id, low, verem, komponensek, elvago_pontok, hidak);
			low[jelen_ind] = min(low[jelen_ind], low[szomszed.ind]);

			if (low[szomszed.ind] > id[jelen_ind]) hidak.push_back(jelen_el_to_szomszed);		// Ha talaltunk hidat

			if ( ((szulo_ind == -1) && gyerekek_szama >= 2) || (szulo_ind != -1 && low[szomszed.ind] >= id[jelen_ind]))		// Ha elvago csomopont
			{
				elvago_pont = true;
				komponensek.resize(komponensek.size() + 1);			// Minden eddig nem latogatott csomopont egy uj komponenst kezd
				while (!verem.empty() && verem.top() != jelen_el_to_szomszed)
				{
					komponensek[komponensek.size() - 1].push_back(verem.top());
					verem.pop();
				}
				komponensek[komponensek.size() - 1].push_back(verem.top());
				verem.pop();
			}
		}
		else if (szomszed.ind != szulo_ind && id[szomszed.ind] < id[jelen_ind])
		{
			low[jelen_ind] = min(low[jelen_ind], id[szomszed.ind]);
			if (id[szomszed.ind] < id[jelen_ind])
			{
				verem.push(jelen_el_to_szomszed);
			}
		}
	}

	// Ket eset van elvago pontnal: 1. Ha nem a jelenlegi bejarasi fa gyokere es 'elvago_pont' teljesult. 2. Ha viszont a bejarasi fa gyokere de van legalabb 2 gyereke
	if (((szulo_ind != -1) && elvago_pont) || ((szulo_ind == -1) && gyerekek_szama >= 2)) elvago_pontok.push_back(jelen_ind);
}

// User functions
bool IranyitatlanGraf::van_kor() const
{
	vector<bool> latogatott(this->n, false), jelenlegi_ut(this->n, false);

	bool letezik_kor = false;
	for (int csomopont = 0; csomopont < this->n && !letezik_kor; ++csomopont)
	{
		if (!latogatott[csomopont])
		{
			letezik_kor = this->melysegi_visit_kor(-1, csomopont, latogatott, jelenlegi_ut);
		}
	}
	return letezik_kor;
}

int IranyitatlanGraf::find_union(int ind, vector<Subset>& subsets) const
{
	if (subsets[ind].parent != ind)
	{
		subsets[ind].parent = this->find_union(subsets[ind].parent, subsets);
	}
	return subsets[ind].parent;
}

void IranyitatlanGraf::union_union(int x, int y, vector<Subset>& subsets) const
{
	int parent_x = this->find_union(x, subsets);
	int parent_y = this->find_union(y, subsets);

	if (subsets[parent_x].rank < subsets[parent_y].rank)
	{
		subsets[parent_x].parent = parent_y;
	}
	else if (subsets[parent_x].rank > subsets[parent_y].rank)
	{
		subsets[parent_y].parent = parent_x;
	}
	else
	{
		subsets[parent_y].parent = parent_x;
		subsets[parent_x].rank++;
	}
}

bool IranyitatlanGraf::van_kor_union_find() const
{
	vector<Subset> subsets;
	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		subsets.push_back(Subset(csomopont, 0));
	}

	for (const El& el : this->el_lista)
	{
		int parent_kezd = find_union(el.kezd, subsets);
		int parent_veg = find_union(el.veg, subsets);

		if (parent_kezd == parent_veg) return true;			// Ha ugyanahoz a set-hez tartoznak akkor talaltunk egy kort

		this->union_union(parent_kezd, parent_veg, subsets);
	}
	return false;
}

pair< list<El>, double > IranyitatlanGraf::Prim_MST() const
{
	/*
		Prim algoritmusa egy MST (minimum spanning tree) algoritmus. Az MST-t es a koltseget teriti.
		Algoritmus mukodese: Kivalasztunk egy kiindulo csomopontot es onnan megjegyezzuk a legrovidebb tavolsagokat
		a tobbi csomopontokba, ha elerhetoek. Minden iteracioban kivalasszuk a legkozelebb levo csomopontot es a hozzavezeto elet.
		Az elso iteracion kivul, minden tobbiben be helyezzuk az adott elet, igy n-1 elunk lesz.
	*/

	int indulas_ind = 0;
	vector<PrimNode> csomopontok(this->n);									// distance = INF, parent = -1, ismeretlen = TRUE
	list<El> minimalis_feszitofa_elek;
	double minimalis_feszitofa = 0;

	csomopontok[indulas_ind].distance = 0;									// Indulasi csomopont (Osszefuggo graf eseten mindegy honnan indulunk)

	for (int elek = 0; elek < this->n; ++elek)								// n - 1 el szukseges, viszont elso iteracioban nem helyez be egy elet sem
	{
		int min_ind = this->get_min_Prim(csomopontok);						// A legkisebb koltsegu csomopontot (A hozza vezeto elet) valassza ki
		csomopontok[min_ind].ismeretlen = false;

		if (min_ind != indulas_ind)
		{
			int parent = csomopontok[min_ind].parent;						// A min_ind parent-je (Ahonnan jovunk)
			minimalis_feszitofa_elek.push_back(El(parent, min_ind, this->szomszedsagi_matrix[min_ind][parent]));	// honnan, hova, suly
			minimalis_feszitofa += this->szomszedsagi_matrix[min_ind][csomopontok[min_ind].parent];
		}

		for (int szomszed_ind = 0; szomszed_ind < this->n; ++szomszed_ind)
		{
			if (csomopontok[szomszed_ind].ismeretlen && this->szomszedsagi_matrix[min_ind][szomszed_ind])				// Ha ismeretlen az a csomopont es letezik ut oda
			{
				if (csomopontok[szomszed_ind].distance > this->szomszedsagi_matrix[min_ind][szomszed_ind])
				{
					csomopontok[szomszed_ind].distance = this->szomszedsagi_matrix[min_ind][szomszed_ind];
					csomopontok[szomszed_ind].parent = min_ind;
				}
			}
		}
	}

	return { minimalis_feszitofa_elek, minimalis_feszitofa };
}

pair< list<El>, double > IranyitatlanGraf::Kruskal_MST() const
{
	/*
		Kruskal algoritmusa egy MST (minimum spanning tree) algoritmus. Az MST-t es a koltseget teriti.
		Az algoritmus mukodese: Kepzeljuk el, hogy eloszor minden csomopont egy onallo komponenshez tartozik, amikor viszont
		egy ellel osszekotunk ket csomopontot, akkor ket onallo komponensbol egy komponens lesz (lenyegtelen melyik veszi fel melyik nevet).
		Eloszor rendezzuk az eleket koltseg szerint novekvo sorrendbe, azutan meg hozzaadjuk a minimalis feszitofahoz az eleket, abban az 
		esetben ha nem hozunk letre egy kort az el hozzaadasaval.
	*/

	list<El> elek(this->el_lista);											// Keszitunk egy masolatot, hogy az eredeti el lista maradjon meg

	elek.sort([](const El& elso, const El& masodik) -> bool { return elso.suly < masodik.suly; });	// Rendezzuk az eleket novekvo sorrendbe

	vector<int> komponensek(this->n);
	for (int i = 0; i < this->n; ++i) komponensek[i] = i;						// Eleinte minden csomopont egy kulon komponenshez tartozik

	// Visszateritendo adatok
	list<El> minimalis_feszitofa;
	double minimalis_feszitofa_koltseg = 0;

	int elek_eddig = 0;
	for(const El& el : elek)												// Atjarjuk az eleket
	{
		if (komponensek[el.kezd] != komponensek[el.veg])
		{
			elek_eddig++;

			minimalis_feszitofa_koltseg += el.suly;
			minimalis_feszitofa.push_back(el);

			int id = komponensek[el.veg];
			for (int k = 0; k < this->n; ++k)									// Frissitjuk a komponenseket
			{
				if (komponensek[k] == id) komponensek[k] = komponensek[el.kezd];
			}

			if (elek_eddig >= this->n - 1) break;								// Ha mar n-1 ele van akkor ez egy fa, tehat megtalaltuk a minimalis feszitofat
		}
	}

	return { minimalis_feszitofa, minimalis_feszitofa_koltseg };
}

pair< list<El>, double > IranyitatlanGraf::Boruvka_MST() const
{
	/*
		Borvuka algoritmusa egy MST (minimum spanning tree) algoritmus. Az MST-t es a koltseget teriti.
		Az algoritmus mukodese: Az algoritmus elso iteracioban megkeresi minden komponenshez (eleinte minden komponens csak egy
		csomopont) tartozo legolcsobb elet. Ezutan pedig egyesiti ezeket az eleket, hogy kevesebb komponenst kapjunk. 
		Ezt a lepest ismeteljuk, amig csak egy nagy komponensunk marad, a minimalis feszitofa.
	*/

	int komponens_db = this->n;								// Eleinte n darab komponensunk (Kulonallo fank van)
	vector<Subset> subsets;
	vector<El> legolcsobb_elek;

	// Vissszateritett minimalis feszitofa es koltsege
	list<El> minimalis_feszitofa;
	double minimalis_feszitofa_koltseg = 0;

	// Inicializaljuk eloszor az n komponenst es legolcsobb eleket
	El kezdeti_el(-1, -1, -1);
	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		subsets.push_back( Subset(csomopont, 0) );
		legolcsobb_elek.push_back(kezdeti_el);
	}

	while (komponens_db > 1)								// Ameddig tobb mint 1 komponensunk van addig osszesitjuk a komponenseket
	{
		for (const El& el : this->el_lista)
		{
			int set1 = this->find_union(el.kezd, subsets);
			int set2 = this->find_union(el.veg, subsets);

			// Ha ugyanahoz a komponenshez tartoznak akkor ignoraljuk
			if (set1 != set2)
			{
				if (legolcsobb_elek[set1] == kezdeti_el || legolcsobb_elek[set1].suly > el.suly)
				{
					legolcsobb_elek[set1] = el;
				}

				if (legolcsobb_elek[set2] == kezdeti_el || legolcsobb_elek[set1].suly > el.suly)
				{
					legolcsobb_elek[set2] = el;
				}
			}
		}

		for (const El& legolcsobb_el : legolcsobb_elek)
		{
			if (legolcsobb_el != kezdeti_el)				// Ha letezik az el
			{
				int set1 = this->find_union(legolcsobb_el.kezd, subsets);
				int set2 = this->find_union(legolcsobb_el.veg, subsets);

				// Ha ugyanahoz a komponenshez tartoznak akkor ignoraljuk
				if (set1 != set2)
				{
					minimalis_feszitofa_koltseg += legolcsobb_el.suly;
					minimalis_feszitofa.push_back(legolcsobb_el);
					this->union_union(set1, set2, subsets);
					komponens_db--;
				}
			}
		}

		// Reset legolcsobb elek
		for (El& el : legolcsobb_elek) el = kezdeti_el;
	}

	return { minimalis_feszitofa, minimalis_feszitofa_koltseg };
}

pair< list<El>, double > IranyitatlanGraf::Forditott_torles_MST() const
{
	/*
		Forditott torles algoritmus egy MST (minimum spanning tree) algoritmus. Az MST-t es a koltseget teriti.
		Az algoritmus mukodese: Az eleket csokkeno sorrendbe rendezzuk es toroljuk azokat az eleket, amely altal a graf
		nem esik szet ket komponensre (Nem hid az adott el).
	*/

	// Vissszateritett minimalis feszitofa es koltsege
	list<El> minimalis_feszitofa;
	double minimalis_feszitofa_koltseg = 0;

									
	IranyitatlanGraf graf(*this);												// Egy el lista masolatot hasznalunk, ugyanis rendeznunk kell az eleket csokkeno sorrendbe
	graf.el_lista.sort([](const El& elso, const El& masodik) -> bool { return elso.suly > masodik.suly; });	// Csokkeno sorrendbe rendezzuk az eleket

	auto it = graf.el_lista.begin();
	while (it != graf.el_lista.end() && graf.el_lista.size() > this->n - 1)		// Ameddig meg lehetseges torolni es tobb mint n-1 elunk van (Egy feszitofa mindig n-1 elbol al)
	{
		list<El> hidak;
		tie(ignore, ignore, hidak) = graf.Tarjan_BiConnect();					// Tarjan Biconnect algoritmusat hasznaljuk ugyanis az egy melysegi bejaras altal megadja nekunk a grafban talalhato hidakat

		// Megnezzuk ha a jelenlegi el hid-e
		// Ha hid akkor nem torolhetjuk, mert a graf szetesik ket komponensre
		bool jelen_el_hid = false;
		for (const El& el : hidak)
		{
			if (el == *it) { jelen_el_hid = true; break; }
		}

		if (jelen_el_hid) it++;													// Ha hid nem toroljuk
		else
		{
			auto it2 = it++;
			graf -= *it2;														// Kivesszuk az adott elet a grafbol
		}
	}

	// Tehat a kapott el lista megegyezik a minimalis feszitofaval
	minimalis_feszitofa = graf.el_lista;
	for (const El& el : minimalis_feszitofa)
	{
		minimalis_feszitofa_koltseg += el.suly;
	}

	return { minimalis_feszitofa, minimalis_feszitofa_koltseg };
}

pair< list<El>, int> IranyitatlanGraf::BFS_SP(int kezd, int veg) const
{
	if (kezd < 0 || veg < 0 || kezd >= this->n || veg >= this->n) throw WrongIndex();

	list<El> ut;															// Visszateritett ut
	int ut_hossz = 0;															// Es hossza

	bool megtalalt = false;
	int jelen_reteg_maradt = 1, kovetkezo_reteg = 0;
	queue<int> sor;
	vector<bool> latogatott(this->n, false);
	vector<int> parent(this->n, -1);

	// Eloszor is behelyezzuk az indulasi koordinatakat
	sor.push(kezd);
	latogatott[kezd] = true;													// Mostmar latogatott

	while (!sor.empty())														// Ameddig meg tudunk keresni es nem talaltuk meg, addig keres
	{
		int jelen_poz = sor.front();											// Jelenlegi pozicio
		sor.pop();

		if (jelen_poz == veg)
		{
			megtalalt = true;
			break;
		}

		for (const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_poz])
		{
			if (!latogatott[szomszed.ind])										// Ha meg nem jartunk azon a szabad helyen
			{
				latogatott[szomszed.ind] = true;
				sor.push(szomszed.ind);
				kovetkezo_reteg++;

				parent[szomszed.ind] = jelen_poz;								// Megjegyezzuk, hogy honnan jottunk
			}
		}

		jelen_reteg_maradt--;
		if (jelen_reteg_maradt == 0)											// Tehat ha elfogyott a jelenlegi 'reteg'-ben a szomszedok, akkor lepunk a kovetkezobe
		{
			jelen_reteg_maradt = kovetkezo_reteg;
			kovetkezo_reteg = 0; ut_hossz++;									// Tehat nott a tavolsag 1-el
		}
	}

	// Utat megkereso resz
	if (megtalalt)
	{
		this->megkeres_utvonal(veg, kezd, ut, parent);

		vector<bool> uton(this->n, false);
		for (const El& el : ut)
		{
			uton[el.kezd] = true;
			uton[el.veg] = true;
		}
	}
	else
	{
		cout << "Nincs megoldas: Nem elerheto a veg pont!";
		ut_hossz = -1;
	}

	return { ut, ut_hossz };
}

tuple< vector< list<El> >, list<int>, list<El> > IranyitatlanGraf::Tarjan_BiConnect() const
{
	/*
		Tarjan BiConnect algoritmusa egy iranyitatlan graf ketszeresen osszefuggo komponenseit, azok szamat, elvago csomopontokat, illetve hidakat hatarozza meg.
		Algoritmus mukodese: Hasonlit a Tarjan StrongConnect algoritmusara, viszont most nem erossen osszefuggo komponenseket kerusunk, hanem ketszeresen
		osszefuggo komponenseket. A hidak es elvago csomopontok megtalalasara kell teljesuljenek bizonyos feltetelek.
	*/

	vector< list<El> > komponensek;							// Ketszeresen osszefuggo komponensek vektor listaja
	list<int> elvago_pontok;									// Elvago pontok listaja
	list<El> hidak;											// Hidak listaja

	int time = 0;
	vector<int> id(this->n, -1), low(this->n, -1);				// id = minden csomopont egy egyedi id-t kap, low = a low-link ertek segit meghatarozni, mely csomopontok tartoznak ugyanahoz a komponenshez
	stack<El> verem;

	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (id[csomopont] == -1)
		{
			this->melysegi_visit_Tarjan_BiConnect(csomopont, -1, time, id, low, verem, komponensek, elvago_pontok, hidak);
		}
	}

	// Ha meg maradt valami a stack-be, akkor az is egy uj ketszeresen osszefuggo komponnens
	if (!verem.empty())
	{
		komponensek.resize(komponensek.size() + 1);
		while (!verem.empty())
		{
			komponensek[komponensek.size() - 1].push_back(verem.top());
			verem.pop();
		}
	}

	return { komponensek, elvago_pontok, hidak };
}

pair< list<El>, double > IranyitatlanGraf::TSP_2_approximate() const
{
	/* MST/DFS / 2_approximate algoritmus implementacio */

	// Meghatarozzuk a grafunk minimalis feszitofajat
	list<El> minimalis_feszitofa;
	tie(minimalis_feszitofa, ignore) = this->Kruskal_MST();

	// Letrehozunk belole egy minimilis feszitofa IranyitatlanGraf objektumot
	IranyitatlanGraf MST(minimalis_feszitofa);

	// Melysegi bejarassal meghatarozzuk a csomopontokat (duplikansok is bekerulnek)
	vector<int> csomopontok = MST.melysegi_bejaras_TSP_2_approximate(1);	// 1-es csomopontbol indulunk

	// Most kitorljuk a duplikansokat (uj vectort hozunk letre)
	vector<bool> volt(this->n, false);
	vector<int> hamilton_kor_csomopontok;
	hamilton_kor_csomopontok.push_back(csomopontok[0]);						// Kezdeti csomopont
	for (vector<int>::const_iterator csp = ++csomopontok.begin(); csp != csomopontok.end(); ++csp)
	{
		if (!volt[*csp])													// A duplikansokat nem helyezzuk bele
		{
			hamilton_kor_csomopontok.push_back(*csp);
			volt[*csp] = true;
		}
	}

	// Letrehozzuk az eleket tartalmazo vector-t es kiszamoljuk a koltseget
	double hamilton_kor_koltseg = 0.0;
	list<El> hamilton_kor;
	for (vector<int>::const_iterator csp = hamilton_kor_csomopontok.begin(); csp != hamilton_kor_csomopontok.end() - 1; ++csp)
	{
		hamilton_kor.push_back(El(*csp, *(csp + 1), szomszedsagi_matrix[*csp][*(csp + 1)]));
		hamilton_kor_koltseg += szomszedsagi_matrix[*csp][*(csp + 1)];
	}

	return { hamilton_kor, hamilton_kor_koltseg };
}
/*--------------------------------------------------------------- Iranyitatlan Graf -------------------------------------------------------------------------*/
