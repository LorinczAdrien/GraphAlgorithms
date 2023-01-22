#include "graf.h"

/*---------------------------------------------------------------- Iranyitott Graf --------------------------------------------------------------------------*/
// Class
IranyitottGraf::IranyitottGraf()
{
	this->n = this->m = 0;
}

IranyitottGraf::IranyitottGraf(string file_name)
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

IranyitottGraf::IranyitottGraf(const vector< vector<double> >& matrix)
{
	if (matrix.size() == 0) throw WrongValue();

	this->n = matrix.size();

	cout << matrix;

	if (matrix.size() == matrix[0].size())						// Ha negyzetes matrix, akkor feltetelezzuk, hogy szomszedsagi matrixot adtak meg
	{
		this->szomszedsagi_matrix = matrix;
		this->m = this->get_elek_szama_from_szomszedsagi_matrix();
		this->szomszedsagi_matrix_to_incidencia_matrix();
		this->incidencia_matrix_to_szomszedsagi_lista();
		this->szomszedsagi_lista_to_el_lista();
	}
	else														// Kulonben incidencia matrix
	{
		this->m = matrix[0].size();
		this->incidencia_matrix = matrix;
		this->incidencia_matrix_to_szomszedsagi_lista();
		this->szomszedsagi_lista_to_el_lista();
		this->el_lista_to_szomszedsagi_matrix();
	}
}

IranyitottGraf::IranyitottGraf(const IranyitottGraf& graf)
{
	this->n = graf.n; this->m = graf.m;

	// El lista masolas
	this->el_lista = graf.el_lista;

	// Szomszedsagi lista masolas
	this->szomszedsagi_matrix.resize(n);
	for (int i = 0; i < this->n; ++i)
	{
		this->szomszedsagi_matrix[i].resize(n);
		for (int j = 0; j < this->n; ++j) this->szomszedsagi_matrix[i][j] = graf.szomszedsagi_matrix[i][j];
	}

	// Incidencia matrix masolas
	this->incidencia_matrix.resize(n);
	for (int i = 0; i < this->n; ++i)
	{
		this->incidencia_matrix[i].resize(m);
		for (int j = 0; j < this->m; ++j) this->incidencia_matrix[i][j] = graf.incidencia_matrix[i][j];
	}

	// Szomszedsagi list masolas
	this->szomszedsagi_lista.resize(n);
	for (int i = 0; i < this->n; ++i)
	{
		this->szomszedsagi_lista[i] = graf.szomszedsagi_lista[i];
	}
}

// Operators =, +, +=, -, -=
IranyitottGraf& IranyitottGraf::operator=(const IranyitottGraf& graf)
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

IranyitottGraf IranyitottGraf::operator +(const El& el)
{
	IranyitottGraf osszeadott(*this);

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

IranyitottGraf IranyitottGraf::operator +(const list<El>& elek)
{
	IranyitottGraf osszeadott(*this);

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
		}
	}

	osszeadott.szomszedsagi_matrix_to_incidencia_matrix();
	osszeadott.incidencia_matrix_to_szomszedsagi_lista();

	return osszeadott;
}

IranyitottGraf IranyitottGraf::operator +(const IranyitottGraf& graf)
{
	IranyitottGraf osszeadott(*this);
	return osszeadott + graf.el_lista;
}

IranyitottGraf& IranyitottGraf::operator +=(const El& el)
{
	return *this = *this + el;
}

IranyitottGraf& IranyitottGraf::operator +=(const list<El>& elek)
{
	return *this = *this + elek;
}

IranyitottGraf& IranyitottGraf::operator +=(const IranyitottGraf& graf)
{
	return *this = *this + graf;
}

IranyitottGraf IranyitottGraf::operator -(const El& el)
{
	IranyitottGraf kivont(*this);

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

IranyitottGraf IranyitottGraf::operator -(const list<El>& elek)
{
	IranyitottGraf kivont(*this);

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

IranyitottGraf IranyitottGraf::operator -(const IranyitottGraf& graf)
{
	IranyitottGraf kivont(*this);
	return kivont - graf.el_lista;
}

IranyitottGraf& IranyitottGraf::operator -=(const El& el)
{
	return *this = *this - el;
}

IranyitottGraf& IranyitottGraf::operator -=(const list<El>& elek)
{
	return *this = *this - elek;
}

IranyitottGraf& IranyitottGraf::operator -=(const IranyitottGraf& graf)
{
	return *this = *this - graf;
}

// Utility
void IranyitottGraf::el_lista_to_szomszedsagi_matrix()
{
	if (this->szomszedsagi_matrix.size() != this->n)
	{
		this->szomszedsagi_matrix.resize(this->n);
		for (int i = 0; i < this->n; ++i) szomszedsagi_matrix[i].resize(this->n);
	}

	for(const El& el : this->el_lista)
	{
		this->szomszedsagi_matrix[el.kezd][el.veg] = el.suly;
	}
}

void IranyitottGraf::el_lista_to_szomszedsagi_lista()
{
	if (this->szomszedsagi_lista.size() != this->n)
	{
		this->szomszedsagi_lista.resize(this->n);							
	}
	else
	{
		this->szomszedsagi_lista = vector< list<SzListaElem> >(n);					// Ha mar letezik, akkor helyebe masolunk egy ures listat
	}

	for (const El& el : this->el_lista)
	{
		this->szomszedsagi_lista[el.kezd].push_back(SzListaElem(el.veg, el.suly));
	}
}

void IranyitottGraf::szomszedsagi_matrix_to_incidencia_matrix()
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
		for (int oszlop_ind = 0; oszlop_ind < this->n; ++oszlop_ind)
		{
			if (this->szomszedsagi_matrix[sor_ind][oszlop_ind])
			{
				this->incidencia_matrix[sor_ind][jelen_el] = this->szomszedsagi_matrix[sor_ind][oszlop_ind];
				this->incidencia_matrix[oszlop_ind][jelen_el] = this->szomszedsagi_matrix[sor_ind][oszlop_ind];

				jelen_el++;															// Csak itt noveljuk, mivel 0-tol indexelunk
			}
		}
	}

}

void IranyitottGraf::incidencia_matrix_to_szomszedsagi_lista()
{
	if (this->szomszedsagi_lista.size() != this->n) this->szomszedsagi_lista.resize(this->n);					// n csomopont
	for (int i = 0; i < this->n; i++) this->szomszedsagi_lista[i].clear();

	for (int oszlop_ind = 0; oszlop_ind < this->m; ++oszlop_ind)
	{
		bool talalt = false;
		int sor_ind = 0, elso_ind = -1, masodik_ind = -1;

		while ((sor_ind < this->n) && !talalt)
		{
			if (this->incidencia_matrix[sor_ind][oszlop_ind] && (elso_ind == -1))
			{
				elso_ind = sor_ind;
			}
			else
			{
				if (this->incidencia_matrix[sor_ind][oszlop_ind] && (masodik_ind == -1))
				{
					masodik_ind = sor_ind;
					talalt = true;								// Tehat megtalaltuk, hogy ez az el melyik ket csomopontot kot ossze (megalhat a while)
				}
			}

			sor_ind++;
		}

		if (this->szomszedsagi_matrix[elso_ind][masodik_ind] && !this->letezik_szomszed( elso_ind, SzListaElem(masodik_ind, this->incidencia_matrix[elso_ind][oszlop_ind]) ) )
		{
			this->szomszedsagi_lista[elso_ind].push_back(SzListaElem(masodik_ind, this->incidencia_matrix[elso_ind][oszlop_ind]));	// Betesszuk a szomszedot es az el sulya
		}
		else if (szomszedsagi_matrix[masodik_ind][elso_ind])
		{
			this->szomszedsagi_lista[masodik_ind].push_back(SzListaElem(elso_ind, this->incidencia_matrix[elso_ind][oszlop_ind]));
		}
	}
}

void IranyitottGraf::szomszedsagi_lista_to_el_lista()
{
	this->el_lista = list<El>();

	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		for(const SzListaElem& szomszed : this->szomszedsagi_lista[sor_ind])
		{
			// Behelyezzuk a csomopont indexet, szomszedjat es az el sulyat
			this->el_lista.push_back( El(sor_ind, szomszed.ind, szomszed.suly));
		}
	}
	this->m = static_cast<int>(this->el_lista.size());
}

void IranyitottGraf::transzponal_el_lista()
{
	for (El& el : this->el_lista)
	{
		swap(el.kezd, el.veg);
	}
}

int IranyitottGraf::get_elek_szama_from_szomszedsagi_matrix() const
{
	int elek_szama = 0;
	for (int sor = 0; sor < this->szomszedsagi_matrix.size(); ++sor)
	{
		for (int oszlop = 0; oszlop < this->szomszedsagi_matrix[sor].size(); ++oszlop)
		{
			if (this->szomszedsagi_matrix[sor][oszlop]) elek_szama++;
		}
	}

	return elek_szama;
}

bool IranyitottGraf::letezik_szomszed(int ind, const SzListaElem& szomszed) const
{
	for (const SzListaElem& ind_szomszed : this->szomszedsagi_lista[ind])
	{
		if (ind_szomszed.ind == szomszed.ind) return true;
	}
	return false;
}

void IranyitottGraf::melysegi_visit_plusz_minusz(int ind, vector<bool>& plusz, vector<bool>& minusz, vector<bool>& latogatott, bool pluszt_keres) const
{
	latogatott[ind] = true;

	// Attol fuggoen, hogy most mivel jelolunk, megjeloljuk a csomopontot
	if (pluszt_keres) plusz[ind] = true;
	else minusz[ind] = true;

	for(const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			this->melysegi_visit_plusz_minusz(szomszed.ind, plusz, minusz, latogatott, pluszt_keres);
		}
	}
}

void IranyitottGraf::melysegi_visit_Kosaraju_verem(int ind, vector<bool>& latogatott, stack<int>& verem) const
{
	latogatott[ind] = true;

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			this->melysegi_visit_Kosaraju_verem(szomszed.ind, latogatott, verem);
		}
	}

	verem.push(ind);												// Miutan minden szomszedjat bejartuk -> a verembe kerul
}

void IranyitottGraf::melysegi_visit_Kosaraju_komponensek(int ind, vector<bool>& latogatott, vector< list<int> >& komponensek) const
{
	latogatott[ind] = true;
	komponensek[komponensek.size() - 1].push_back(ind);				// Behelyezzuk a csomopontot a megfelelo komponensbe

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			this->melysegi_visit_Kosaraju_komponensek(szomszed.ind, latogatott, komponensek);
		}
	}
}

void IranyitottGraf::melysegi_visit_Tarjan(int jelen_ind, int& time, vector< list<int> >& komponensek, stack<int>& verem, vector<int>& id, vector<int>& low, vector<bool>& vermen) const
{
	id[jelen_ind] = low[jelen_ind] = ++time;											// Belepeskor a low es id time megegyezik
	verem.push(jelen_ind);
	vermen[jelen_ind] = true;															// Jelenleg a vermen van a csomopont			

	for(const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_ind])
	{
		if (id[szomszed.ind] == -1)														// Ha meg nem voltunk ott
		{
			this->melysegi_visit_Tarjan(szomszed.ind, time, komponensek, verem, id, low, vermen);
			low[jelen_ind] = min(low[jelen_ind], low[szomszed.ind]);					// Minimumot veszi fel a visszateresnel
		}
		else if ((id[szomszed.ind] < id[jelen_ind]) && (vermen[szomszed.ind]))
		{
			low[jelen_ind] = min(low[jelen_ind], id[szomszed.ind]);
		}
	}

	if (low[jelen_ind] == id[jelen_ind])
	{
		komponensek.resize(komponensek.size() + 1);										// Tehat talaltunk egy uj osszefuggo komponenst
		while ((!verem.empty()) && (id[verem.top()] >= id[jelen_ind]))					// Az ugyanahhoz a komponenshez tartozo csucspontokat eltavolitjuk a stack-rol
		{
			komponensek[komponensek.size() - 1].push_back(verem.top());
			vermen[verem.top()] = false;
			verem.pop();
		}
	}
}

void IranyitottGraf::init_Floyd_Warshall(vector< vector<double> >& tavolsag, vector< vector<int> >& parent) const
{
	for (int sor = 0; sor < this->n; ++sor)
	{
		for (int oszlop = 0; oszlop < this->n; ++oszlop)
		{
			if (this->szomszedsagi_matrix[sor][oszlop] != 0)							// Ha nem vegtelen a tavolsag bemasoljuk
			{
				tavolsag[sor][oszlop] = this->szomszedsagi_matrix[sor][oszlop];
				parent[sor][oszlop] = oszlop;
			}
		}
	}

	for (int atlo = 0; atlo < this->n; ++atlo) tavolsag[atlo][atlo] = 0;				// Foatlot lenullazuk
}

vector<int> IranyitottGraf::szelessegi_bejaras_Edmonds_csomopontok(int indulas_ind, vector<vector<double> >& rezidualis_graf)
{
	if (indulas_ind < 0 || indulas_ind > this->n) throw WrongIndex();

	vector<bool> latogatott(this->n, false);
	vector<int> erintett_csomopontok;
	queue<int> sor;

	erintett_csomopontok.push_back(indulas_ind);											// alapbol benne van ahonnan indulunk
	sor.push(indulas_ind);
	latogatott[indulas_ind] = true;

	while (!sor.empty())
	{
		int jelen_ind = sor.front();
		sor.pop();

		for (int csomopont_ind = 0; csomopont_ind < this->n; ++csomopont_ind)
		{
			if (!latogatott[csomopont_ind] && (rezidualis_graf[jelen_ind][csomopont_ind]))
			{
				erintett_csomopontok.push_back(csomopont_ind);								// tehat eltudtuk erni
				sor.push(csomopont_ind);
				latogatott[csomopont_ind] = true;
			}
		}
	}

	return erintett_csomopontok;
}

bool IranyitottGraf::szelessegi_bejaras_Edmonds(int forras, int nyelo, vector<int>& parent, vector<vector<double>>& rezidualis_graf)
{
	if (forras < 0 || forras > this->n) throw WrongIndex();

	queue<int> sor;
	vector<bool> latogatott(this->n, false);

	sor.push(forras);
	latogatott[forras] = true;
	parent[forras] = -1;

	while (!sor.empty())
	{
		int jelen_ind = sor.front();
		sor.pop();

		for (int csomopont_ind = 0; csomopont_ind < this->n; ++csomopont_ind)
		{
			if (!latogatott[csomopont_ind] && (rezidualis_graf[jelen_ind][csomopont_ind]))	// Ha nem latogatt meg es letezik ut
			{
				if (csomopont_ind == nyelo)													// Ha elertunk a nyelohoz akkor vissza terhetunk
				{
					parent[csomopont_ind] = jelen_ind;
					return true;
				}

				sor.push(csomopont_ind);
				latogatott[csomopont_ind] = true;
				parent[csomopont_ind] = jelen_ind;
			}
		}
	}

	return false;
}

void IranyitottGraf::Pumpalo_init(int forras, vector<PumpaloNode>& csomopontok, vector<vector<double>>& rezidualis_graf)
{
	csomopontok[forras].magassag = this->n;													// Forras magassaga: n

	int szomszed_darab = 0;
	for(const SzListaElem& szomszed : this->szomszedsagi_lista[forras])
	{
		double kapacitas = this->szomszedsagi_matrix[forras][szomszed_darab++];

		// kipumpaljuk amit tudunk a forrasbol
		csomopontok[szomszed.ind].tobbletfolyam = kapacitas;
		rezidualis_graf[forras][szomszed.ind] -= kapacitas;
		rezidualis_graf[szomszed.ind][forras] += kapacitas;
	}
}

void IranyitottGraf::Pumpalas(const El& el, vector<PumpaloNode>& csomopontok, vector<vector<double>>& rezidualis_graf)
{
	// kiszamitjuk mennyit tudunk pumpalni az adott elen
	double pumpalhato = min(rezidualis_graf[el.kezd][el.veg], csomopontok[el.kezd].tobbletfolyam);

	// a rezidualis grafot frisstijuk
	rezidualis_graf[el.kezd][el.veg] -= pumpalhato;
	rezidualis_graf[el.veg][el.kezd] += pumpalhato;

	// a csomopontok tobbletfolyamat frissitjuk
	csomopontok[el.kezd].tobbletfolyam -= pumpalhato;
	csomopontok[el.veg].tobbletfolyam += pumpalhato;
}

void IranyitottGraf::Emeles(int emelt, vector<PumpaloNode>& csomopontok, vector<vector<double>>& rezidualis_graf)
{
	int m = INT32_MAX;
	for (int csomopont_ind = 0; csomopont_ind < this->n; ++csomopont_ind)
	{
		if (rezidualis_graf[emelt][csomopont_ind])										// Ha letezik el
		{
			m = min(m, csomopontok[csomopont_ind].magassag);
		}
	}

	csomopontok[emelt].magassag = m + 1;
}

bool IranyitottGraf::pumpalhato_el(const El& el, int forras, int nyelo, const vector<PumpaloNode>& csomopontok, const vector<vector<double>>& rezidualis_graf)
{
	if ((el.kezd != forras) && (el.kezd != nyelo) && (csomopontok[el.kezd].magassag == csomopontok[el.veg].magassag + 1))		// Ha nem forras, nem nyelo es ha megengedett el
	{
		if (csomopontok[el.kezd].tobbletfolyam > 0 && rezidualis_graf[el.kezd][el.veg] > 0)										// Kezd aktiv csomopont es van amit pumpalni
		{
			return true;
		}
	}
	return false;
}

bool IranyitottGraf::emelheto_csomopont(int csomopont, int forras, int nyelo, const vector<PumpaloNode>& csomopontok, const vector<vector<double>>& rezidualis_graf)
{
	if (csomopontok[csomopont].magassag < this->n)
	{
		if ((csomopont != forras) && (csomopont != nyelo) && (csomopontok[csomopont].tobbletfolyam > 0))
		{
			for (int szomszed = 0; szomszed < this->n; ++szomszed)
			{
				if ((szomszed != forras) && (szomszed != nyelo) && (rezidualis_graf[csomopont][szomszed] > 0))
				{
					if (csomopontok[szomszed].magassag < csomopontok[csomopont].magassag)
					{
						return false;
					}
				}
			}
			return true;
		}
	}
	return false;
}

// User functions
vector< vector<bool> > IranyitottGraf::tranzitiv_lezaras_matrix() const
{
	vector < vector<bool> > letezik_ut(this->n, vector<bool>(this->n, false));

	// Matrix inicializalas
	for (int sor = 0; sor < this->n; ++sor)
	{
		for (int oszlop = 0; oszlop < this->n; ++oszlop)
		{
			if (this->szomszedsagi_matrix[sor][oszlop] != 0)							// Ha nem vegtelen a tavolsag bemasoljuk
			{
				letezik_ut[sor][oszlop] = true;
			}
		}
	}
	for (int atlo = 0; atlo < this->n; ++atlo) letezik_ut[atlo][atlo] = true;			// Foatlo kulon eset

	for (int k = 0; k < this->n; ++k)
	{
		for (int sor = 0; sor < this->n; ++sor)
		{
			for (int oszlop = 0; oszlop < this->n; ++oszlop)
			{
				// Ha letezik ut sor->k es k->oszlop, akkor letezik a sor->oszlop ut is
				letezik_ut[sor][oszlop] = letezik_ut[sor][oszlop] || (letezik_ut[sor][k] && letezik_ut[k][oszlop]);
			}
		}
	}

	return letezik_ut;
}

IranyitottGraf IranyitottGraf::tranzitiv_lezaras() const
{
	vector< vector<double> > tavolsag_matrix;
	tie(tavolsag_matrix, ignore) = this->Floyd_Warshall();

	for (int i = 0; i < tavolsag_matrix.size(); ++i)
	{
		for (int j = 0; j < tavolsag_matrix[i].size(); ++j)
		{
			if (tavolsag_matrix[i][j] == DBL_MAX / 2) tavolsag_matrix[i][j] = 0;	// Floyd-Warshall-nal a vegtelen jelzi, hogy nincs ut, viszont a szomszedsagi matrixban mi maskepp taroltuk
		}
	}

	return IranyitottGraf(tavolsag_matrix);
}

vector< list<int> > IranyitottGraf::plusz_minusz_algoritmus() const
{
	/*
		A 'plusz-minusz' algoritmus egy iranyitott graf teljesen osszefuggo komponenseit teriti. (Tartalmazo matrixot)
		Az algoritmus mukodese: Minden eddig nem egy komponenshez tartozo csomopontbol inditunk egy melysegi bejarast
		ahol eloszor is megjeloljuk az elerheto csomopontokat 'plusz'-al. Ezutan transzponaljuk az el listat (Last transzponal_el_lista() ).
		Transzponalas utan pedig ugyancsak bejarjuk a grafot az adott csomopontbol es megjeloljuk az elerheto csomopontokat 'minusz'-al.
		Vegul pedig minden olyan csomopont amely 'plusz-minusz' es nem tartozik meg egy komponenshez sem azokat egy komponensbe soroljuk.
		Ismeteljuk az algoritmust minden nem komponensben levo csomopontra.
	*/

	IranyitottGraf graf(*this);											// Letrehozunk egy ideiglenes objektumot, hogy hasznalhassuk a 'const'
																		// minositot a fuggvenyen (A transzponal el lista nem 'const' fuggveny)

	vector< list<int> > komponensek;
	vector<bool> plusz(graf.n, false), minusz(graf.n, false), latogatott(graf.n, false), komponensben_van(graf.n, false);

	for(int csomopont = 0; csomopont < graf.n; ++csomopont)
	{
		if (!komponensben_van[csomopont])
		{
			latogatott = plusz = minusz = vector<bool>(graf.n, false);

			// Bejarjuk a grafot az adott csomopontbol eloszor 'plusz'-al jelolve minden olyan csomopontot ahova eljutunk
			graf.melysegi_visit_plusz_minusz(csomopont, plusz, minusz, latogatott, true);

			graf.transzponal_el_lista();								// Transzponaljuk az el listat (Megforditjuk)
			graf.el_lista_to_szomszedsagi_lista();						// Frissit szomszedsagi lista

			latogatott = vector<bool>(graf.n, false);

			// Bejarjuk a grafot az adott csomopontbol 'minusz'-al jelolve minden olyan csomopontot ahova eljutunk
			graf.melysegi_visit_plusz_minusz(csomopont, plusz, minusz, latogatott, false);

			graf.transzponal_el_lista();								// Helyrehozzuk az el listat
			graf.el_lista_to_szomszedsagi_lista();						// Frissit szomszedsagi lista

			bool van_uj_komponens = false;
			for (int i = 0; i < graf.n; ++i)
			{
				if (plusz[i] && minusz[i] && !komponensben_van[i])		// Ha +- es ha meg nem tartozik egy komponenshez sem
				{
					if (!van_uj_komponens) komponensek.resize(komponensek.size() + 1);

					van_uj_komponens = true;
					komponensben_van[i] = true;							// Tehat az i-dik csomopont mostmar egy komponensben van

					komponensek[komponensek.size() - 1].push_back(i);	// A csomopont bekerul a megfelelo komponensbe
				}
			}
		}
	}

	return komponensek;
}

vector< list<int> > IranyitottGraf::Kosaraju_algoritmusa() const
{
	/*
		Kosaraju algoritmusa egy iranyitott graf teljesen osszefuggo komponenseit teriti. (Tartalmazo matrixot)
		Az algoritmus mukodese: Eloszor bejarjuk a grafot es a csomopontokat egy verembe helyezzuk. (Amikor bejartuk
		a csomopont osszes szomszedjat akkor kerul a verembe)
		Ezutan bejarjuk a transzponalt grafot a verem szerint.
	*/

	IranyitottGraf graf(*this);											// Letrehozunk egy ideiglenes objektumot, hogy hasznalhassuk a 'const'
																		// minositot a fuggvenyen (A transzponal el lista nem 'const' fuggveny)

	vector< list<int> > komponensek;									// Visszateritett komponens matrix
	
	stack<int> verem;
	vector<bool> latogatott(graf.n, false);

	// 1. Bejarjuk a vermet es a csomopontokat egy verembe helyezzuk (Mikor minden szomszedjat bejartuk akkor kerul a verembe)
	for (int csomopont = 0; csomopont < graf.n; ++csomopont)
	{
		if (!latogatott[csomopont])
		{
			graf.melysegi_visit_Kosaraju_verem(csomopont, latogatott, verem);
		}
	}

	// 2. Bejarjuk a ! transzponalt ! grafot a verem szerint
	graf.transzponal_el_lista();
	graf.el_lista_to_szomszedsagi_lista();

	latogatott = vector<bool>(graf.n, false);
	while ( !verem.empty() )
	{
		int jelen_csomopont = verem.top();
		verem.pop();

		if (!latogatott[jelen_csomopont])
		{
			komponensek.resize(komponensek.size() + 1);					// Tehat talaltunk egy uj komponenst
			graf.melysegi_visit_Kosaraju_komponensek(jelen_csomopont, latogatott, komponensek);
		}
	}

	return komponensek;
}

vector< list<int> > IranyitottGraf::Tarjan_StrongConnect() const
{
	/*
		Tarjan StrongConnect algoritmusa egy iranyitott graf teljesen osszefuggo komponenseit teriti. (Tartalmazo matrixot)
		Az algoritmus mukodese: Egy DFS bejaras alatt minden csomoponthoz egy egyedi 'id'-t rendelunk, visszateresnel, illetve mar latogatott
		csomopontnal minimumot szamolunk (Lasd rekurziv hivas). 1. eset: Ha a szomszedhoz meg nem rendeltunk egy 'id'-t akkor elvegezzuk a rekurziv
		hivast a szomszedon es azutan 'low[jelen_ind] = min(low[jelen_ind], low[szomszed.ind])'. 2. eset: Ha mar ismert a csomopont akkor viszont csak
		akkor szamitunk minimumot, ha a szomszed 'id'-ja kisebb es a vermen van. 'low[jelen_ind] = min(low[jelen_ind], id[szomszed.ind])'.
	*/

	int time = 0;
	vector< list<int> > komponensek;					// A komponenseket tarolo lista vektor (Minden sor egy komponens)
	vector<int> id(this->n, -1), low(this->n, -1);		// id = tulajdonkeppen minden csompont egy egyedi id-t kap, low = a low-link ertek segit meghatarozni, mely csomopontok tartoznak ugyanahoz a komponenshez
	vector<bool> vermen(this->n, false);				// Egy vermet es egy verem bool vektort is hasznalunk, hogy O(1) idoben le tudjuk kerni, ha egy csomopont a vermen van
	stack<int> verem;

	// Meglatogatunk minden csomopontot es szomszedait amelyek meg nem kaptak egyedi 'id' value-t
	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (id[csomopont] == -1)
		{
			this->melysegi_visit_Tarjan(csomopont, time, komponensek, verem, id, low, vermen);
		}
	}

	return komponensek;
}

pair< vector< vector<double> >, vector< vector<int> > > IranyitottGraf::Johnson()
{
	// 1. lepes: 'hozzaadunk' egy fiktiv csucsot a grafhoz
	int volt_n = this->n;														// Lement csomopontok szama
	this->n++;
	int fiktiv_csomopont = this->n - 1;

	// 2. lepes: A fiktiv csucsot osszekotjuk a tobbi csuccsal 0-s koltsegu ellel
	int volt_m = this->m;														// Lement elek szama
	list<El> lementett_el_lista(this->el_lista);								// Lementjuk az eredeti listat
	for (int csomopont = 0; csomopont < volt_n; ++csomopont) this->el_lista.push_back(El(fiktiv_csomopont, csomopont, 0));		// n darab 0 koltsegu el ami fiktivbol megy egy csomopontba
	this->m += volt_n;

	// 3. lepes: Lefuttatjuk a Bellman-Ford algoritmust 'fiktiv_csomopont' indulasi index-el
	vector<double> Bellman_tavolsag;
	tie(Bellman_tavolsag, ignore) = this->Bellman_Ford_SP(fiktiv_csomopont + 1);// +1 mert alapbol kivon egyet (es ignore, jelzi a tie-nak, hogy a parent vector-t torolje)

	if (Bellman_tavolsag.size() == 0)											// Ha van benne negativ kor akkor megalitjuk az algoritmust
	{
		// Visszaallitjuk az eredeti ertekeket
		this->n = volt_n;
		this->m = volt_m;
		this->el_lista = lementett_el_lista;
		throw GraphHasNegativeCycle();
	}
	else
	{
		// 4. lepes: Minden elet ujrasulyozunk
		for (El& el : this->el_lista)
		{
			el.suly = el.suly + Bellman_tavolsag[el.kezd] - Bellman_tavolsag[el.veg];
		}

		// 5. lepes: Toroljuk a fiktiv csomopontot es futtatjuk Dijkstra algoritmusat minden csomopontra
		this->n = volt_n;
		this->m = volt_m;
		this->el_lista.resize(this->m);											// Toroljuk az eddig hozza adott eleket

		vector< vector<double> > tavolsag(this->n, vector<double>(this->n, DBL_MAX));
		vector< vector<int> > parent(this->n, vector<int>(this->n, -1));

		for (int csomopont = 0; csomopont < this->n; ++csomopont)
		{
			tie(tavolsag[csomopont], parent[csomopont]) = this->Dijkstra_SP(csomopont + 1);
		}

		// Visszaallitjuk az eredeti el listat
		this->el_lista = lementett_el_lista;

		return { tavolsag, parent };
	}
}

pair< vector< vector<double> >, vector< vector<int> > > IranyitottGraf::Floyd_Warshall() const
{
	vector<vector<double>> tavolsag(this->n, vector<double>(this->n, DBL_MAX / 2));				// Megprobaljuk elkerulni az overflow lehetoseget ezert DBL_MAX/2 -t hasznalunk
	vector<vector<int>> parent(this->n, vector<int>(this->n, -1));

	this->init_Floyd_Warshall(tavolsag, parent);												// Inicializaljuk a tavolsag es parent matrixot

	for (int k = 0; k < this->n; ++k)
	{
		for (int sor = 0; sor < this->n; ++sor)
		{
			for (int oszlop = 0; oszlop < this->n; ++oszlop)
			{
				if (tavolsag[sor][k] + tavolsag[k][oszlop] < tavolsag[sor][oszlop])				// Ha talaltunk egy jobb utat adott k-n keresztul, akkor frissitjuk
				{
					tavolsag[sor][oszlop] = tavolsag[sor][k] + tavolsag[k][oszlop];
					parent[sor][oszlop] = parent[sor][k];
				}
			}
		}
	}

	return { tavolsag, parent };
}

tuple<double, vector<int>, vector<int> > IranyitottGraf::Edmonds_Karp(int forras, int nyelo)
{
	forras--; nyelo--;																			// 0-tol indexeles
	vector<vector<double>> rezidualis_graf = this->szomszedsagi_matrix;							// Letrehozunk egy rezidualis grafot, ami eloszor az eredeti szomszedsagi matrix ertekeit tartalmazza
	vector<int> parent(this->n, -1);
	vector<int> min_vagat_halmaz1, min_vagat_halmaz2;

	double maximalis_folyam = 0;

	while (this->szelessegi_bejaras_Edmonds(forras, nyelo, parent, rezidualis_graf))			// Amig eltudunk jutni a nyelohoz a rezidualis grafban
	{
		double jelen_maximalis_folyam = DBL_MAX;

		for (int ind = nyelo; ind != forras; ind = parent[ind])
		{
			int szulo = parent[ind];
			jelen_maximalis_folyam = min(jelen_maximalis_folyam, rezidualis_graf[szulo][ind]);
		}

		for (int ind = nyelo; ind != forras; ind = parent[ind])
		{
			int szulo = parent[ind];
			rezidualis_graf[szulo][ind] -= jelen_maximalis_folyam;
			rezidualis_graf[ind][szulo] += jelen_maximalis_folyam;
		}

		maximalis_folyam += jelen_maximalis_folyam;
	}

	min_vagat_halmaz1 = this->szelessegi_bejaras_Edmonds_csomopontok(forras, rezidualis_graf);	// Minimalis vagat elso halmaza (Amelyek elerhetoek a forrasbol a rezidualis graf altal)

	vector<bool> eleme_halmaz(this->n, false);
	for (int csomopont = 0; csomopont < min_vagat_halmaz1.size(); ++csomopont)					// Megnezzuk melyik csomopontokat hasznaltunk mar
	{
		eleme_halmaz[min_vagat_halmaz1[csomopont]] = true;
	}

	for (int csomopont = 0; csomopont < this->n; ++csomopont)									// Minimalis vagat masodik halmazasnak meghatarozasa (amelyek elerhetoek a nyelobol)
	{
		if (!eleme_halmaz[csomopont]) min_vagat_halmaz2.push_back(csomopont);
	}

	return { maximalis_folyam, min_vagat_halmaz1, min_vagat_halmaz2 };
}

double IranyitottGraf::Pumpalo_algoritmus(int forras, int nyelo)
{
	forras--; nyelo--;																						// 0-tol indexelunk
	vector<PumpaloNode> csomopontok(this->n);																// Alapertelmezetten: magassag = 0, ertek = 0.0
	vector<vector<double>> rezidualis_graf(this->szomszedsagi_matrix);										// Eleinte a rezidualis graf megegyezik a szomszedsagi_matrix-al
	this->Pumpalo_init(forras, csomopontok, rezidualis_graf);

	bool volt_muvelet = true;
	while (volt_muvelet)
	{
		volt_muvelet = false;

		for (const El& el : this->el_lista)																// Megnezzuk ha valamelyik elen keresztul tudunk e pumpalni
		{
			if (this->pumpalhato_el(el, forras, nyelo, csomopontok, rezidualis_graf))
			{
				this->Pumpalas(el, csomopontok, rezidualis_graf);
				volt_muvelet = true;
				break;
			}

			if (this->pumpalhato_el(El(el.veg, el.kezd), forras, nyelo, csomopontok, rezidualis_graf))	// Megnezzuk a rezidualis eleket is
			{
				Pumpalas({ el.veg, el.kezd }, csomopontok, rezidualis_graf);
				volt_muvelet = true;
				break;
			}
		}

		if (volt_muvelet) continue;																			// Ha sikerult pumpalni akkor eloszor megprobalunk megint pumpalni

		for (int csomopont = 0; csomopont < this->n; ++csomopont)
		{
			if (this->emelheto_csomopont(csomopont, forras, nyelo, csomopontok, rezidualis_graf))			// Minden csomopontra megnezzuk ha tudunk-e emelni rajta
			{
				this->Emeles(csomopont, csomopontok, rezidualis_graf);
				volt_muvelet = true;
				break;
			}
		}
	}

	return csomopontok[nyelo].tobbletfolyam;														// a nyelo tobbletfolyama tartalmazza a maximalis folyam erteket
}

vector<TevekenysegNode> IranyitottGraf::Kritikus_ut_masodik_modell(const vector<double>& vegrehajtasi_idok)
{
	vector<TevekenysegNode> csomopontok;
	for (const double& ido : vegrehajtasi_idok) csomopontok.push_back(TevekenysegNode(ido));

	vector<int> topologikus_sorrend = this->topologiai_rendezes();								// Korabbi hazikban megirt topologiai rendezest felhasznaljuk

	// Inicializaljuk a elso csomopontot (topologikusan elso) legkorabbi ertekeit
	int topo_elso = topologikus_sorrend[0];
	csomopontok[topo_elso].legkorabbi_kezdes = 0;
	csomopontok[topo_elso].legkorabbi_befejezodes = csomopontok[topo_elso].vegrehajtasi_ido;

	// Bejarjuk topologikus sorrendben elolrol a grafot es kiszamitjuk a legkorabbi kezdes/befejezodes ertekeit
	for (int csomopont = 1; csomopont < this->n; ++csomopont)
	{
		int jelen_ind = topologikus_sorrend[csomopont];

		for (int bemeno_szomszed = 0; bemeno_szomszed < this->n; ++bemeno_szomszed)				// Megnezzuk minden csomopont eseten, ha 'bemegy' az adott csomopontba
		{
			if (this->szomszedsagi_matrix[bemeno_szomszed][jelen_ind])							// Ha letezik bemeno ut a szomszedbol a 'jelen_ind'-be
			{
				csomopontok[jelen_ind].legkorabbi_kezdes = max(csomopontok[jelen_ind].legkorabbi_kezdes, csomopontok[bemeno_szomszed].legkorabbi_befejezodes);
			}
		}

		csomopontok[jelen_ind].legkorabbi_befejezodes = csomopontok[jelen_ind].legkorabbi_kezdes + csomopontok[jelen_ind].vegrehajtasi_ido;
	}

	// Inicializaljuk az utolso csomopont (topologikusan utolso) legkesobbi ertekeit
	int topo_utolso = topologikus_sorrend[n - 1];
	csomopontok[topo_utolso].legkesobbi_befejezodes = csomopontok[topo_utolso].legkorabbi_befejezodes;
	csomopontok[topo_utolso].legkesobbi_kezdes = csomopontok[topo_utolso].legkesobbi_befejezodes - csomopontok[topo_utolso].vegrehajtasi_ido;

	// Bejarjuk topologikus sorrendben hatulrol a grafot es kiszamitjuk a legkesobbi kezdes/befejezodes ertekeit
	for (int csomopont = this->n - 2; csomopont >= 0; --csomopont)
	{
		int jelen_ind = topologikus_sorrend[csomopont];

		for (const SzListaElem& kimeno_szomszed : this->szomszedsagi_lista[jelen_ind])	// Hasznalhatjuk a szomszedsagi listat ugyanis most a 'jelen_ind' csucsbol kimeno eleket nezzuk es azon csomopontjait
		{
			csomopontok[jelen_ind].legkesobbi_befejezodes = min(csomopontok[jelen_ind].legkesobbi_befejezodes, csomopontok[kimeno_szomszed.ind].legkesobbi_kezdes);
		}

		csomopontok[jelen_ind].legkesobbi_kezdes = csomopontok[jelen_ind].legkesobbi_befejezodes - csomopontok[jelen_ind].vegrehajtasi_ido;
	}

	return csomopontok;
}
/*---------------------------------------------------------------- Iranyitott Graf --------------------------------------------------------------------------*/