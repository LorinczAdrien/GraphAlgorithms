#include "graf.h"

/*-------------------------------------------------------------------- Graf ---------------------------------------------------------------------------------*/
// Utility
void Graf::melysegi_visit(int ind, vector<bool>& latogatott, vector<int>& bejart) const
{
	latogatott[ind] = true;

	for(const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			bejart.push_back(szomszed.ind);
			this->melysegi_visit(szomszed.ind, latogatott, bejart);
		}
	}
}

bool Graf::melysegi_visit_kor(int ind, vector<bool>& latogatott, vector<bool>& jelenlegi_ut) const
{
	bool van_kor = false;
	latogatott[ind] = true;
	jelenlegi_ut[ind] = true;

	for(const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (jelenlegi_ut[szomszed.ind])						// Ha a jelenlegi utunkban mar szerepelt ez a csomopont (jelenleg egy szomszed), akkor talaltunk egy kort
		{
			van_kor = true;
			break;
		}
		else
		{
			if (!latogatott[szomszed.ind])
			{
				van_kor = this->melysegi_visit_kor(szomszed.ind, latogatott, jelenlegi_ut);
			}
		}
	}

	jelenlegi_ut[ind] = false;								// Visszalepeskor frissitjuk a jelenlegi utat (mint backtrack eseten)

	return van_kor;
}

void Graf::nodes_to_node_recursive(const vector<int>& parent, vector<int>& ut, int ind) const
{
	if (parent[ind] != -1)
	{
		this->nodes_to_node_recursive(parent, ut, parent[ind]);
		ut.push_back(ind);
	}
	else
	{
		ut.push_back(ind);
	}
}

void Graf::topologiai_visit(int ind, vector<bool>& latogatott, vector<int>& topo) const
{
	latogatott[ind] = true;

	for (const SzListaElem& szomszed : this->szomszedsagi_lista[ind])
	{
		if (!latogatott[szomszed.ind])
		{
			this->topologiai_visit(szomszed.ind, latogatott, topo);
		}
	}

	topo.push_back(ind);
}

vector<int> Graf::find_path_of_nodes_to_node(const vector<int>& parent, int ind) const
{
	vector<int> ut;

	if (parent[ind] != -1)														// Ha nem letezik ut -> ures vector-t teritunk
	{
		this->nodes_to_node_recursive(parent, ut, ind);
	}

	return ut;
}

// User functions
void Graf::torol_csomopont(int ind)
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

	this->el_lista_to_szomszedsagi_matrix();
	this->szomszedsagi_matrix_to_incidencia_matrix();
	this->incidencia_matrix_to_szomszedsagi_lista();
}

void Graf::kiir_el_lista() const
{
	cout << "El lista: " << '\n';
	for(const El& el : this->el_lista)
	{
		cout << el.kezd + 1 << ' ' << el.veg + 1 << ' ' << el.suly << '\n';
	}
	cout << '\n';
}

void Graf::kiir_szomszedsagi_matrix() const
{
	cout << "Szomszedsagi matrix: " << '\n';
	for (int i = 0; i < this->n; ++i)
	{
		cout << i + 1 << ": ";
		for (int j = 0; j < this->n; ++j)
		{
			cout << this->szomszedsagi_matrix[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
}

void Graf::kiir_incidencia_matrix() const
{
	cout << "Incidencia matrix: " << endl;
	for (int i = 0; i < this->n; ++i)
	{
		cout << i + 1 << ": ";
		for (int j = 0; j < this->m; ++j) cout << this->incidencia_matrix[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';
}

void Graf::kiir_szomszedsagi_lista() const
{
	cout << "Szomszedsagi lista: " << '\n';
	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		cout << sor_ind + 1 << ": ";
		for(const SzListaElem& szomszed : this->szomszedsagi_lista[sor_ind])
		{
			cout << szomszed.ind + 1 << '(' << szomszed.suly << ") ";
		}
		cout << '\n';
	}
	cout << '\n';
}

vector<int> Graf::szelessegi_bejaras(int indulas_ind) const
{
	if (indulas_ind < 0 || indulas_ind > this->n) throw WrongIndex();

	/*
		A szélességi bejárás egy gráf bejárási lehetõség.
		A lényege, hogy a csomópontokat egy 'hullám' féle képpen járja be. Ez azt jelenti, hogy a közelebb levõ csomópontokat majdnem egyszerre találja meg.
		Egy prioritás sort használ amely által garantálja, hogy elõször az 1 távolságra levõ csomópontokat találja meg azután 2, 3, 4... .
		Ezért használható súlyozatlan gráfokban legrövidebb út megtalálására.
		Egy std::vector<int> típusú tárolót térit, amely tartalmazza a bejárt sorrendben a csomópontokat.
	*/

	vector<int> bejart;
	vector<bool> latogatott(this->n, false);
	queue<int> sor;

	if (indulas_ind)							// Ha nem az alapertelmezett eset (Tehat megadtunk egy indexet)
	{
		indulas_ind--;
		sor.push(indulas_ind);
		latogatott[indulas_ind] = true;
		indulas_ind++;
	}

	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (indulas_ind == 0 && !latogatott[csomopont])
		{
			sor.push(csomopont);
			latogatott[csomopont] = true;
		}

		while (!sor.empty())
		{
			int jelen_ind = sor.front();
			sor.pop();

			bejart.push_back(jelen_ind);

			for(const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_ind])
			{
				if (!latogatott[szomszed.ind])
				{
					latogatott[szomszed.ind] = true;
					sor.push(szomszed.ind);
				}
			}
		}

		if (indulas_ind) return bejart;		// Alapertelmezett esetben nem lep ki
	}

	return bejart;
}

vector<int> Graf::melysegi_bejaras(int indulas_ind) const
{
	if (indulas_ind < 0 || indulas_ind > this->n) throw WrongIndex();

	/*
		A mélységi bejárás egy gráf bejárási lehetõség.
		A lényege az, hogy a BFS-el ellentétben elõször a kiíndulasi csomópont egy szomszédját válassza ki, ezután pedig ennek a szomszédnak az egyik szomszédját es így tovább.
		Ezáltal egy DFS hivás elõször mélyre halad és ezeket az utakat találja meg.
		Egy std::vector<int> típusú tárolót térit, amely tartalmazza a bejárt sorrendben a csomópontokat.
	*/

	vector<int> bejart;
	vector<bool> latogatott(this->n, false);

	if (indulas_ind)										// Ha megadtunk egy indexet, akkor csak abbol az indexbol jarjuk be a grafot
	{
		indulas_ind--;
		bejart.push_back(indulas_ind);						// Indulasi csomopontot berakjuk

		this->melysegi_visit(indulas_ind, latogatott, bejart);
	}
	else													// Ellenkezo esetben pedig minden csomopontbol bejarjuk a grafot
	{
		for (int csomopont = 0; csomopont < this->n; ++csomopont)
		{
			if (!latogatott[csomopont])
			{
				bejart.push_back(csomopont);
				this->melysegi_visit(csomopont, latogatott, bejart);
			}
		}
	}

	return bejart;
}

bool Graf::regularis() const
{
	int r_ertek = 0;												// Az r erteket az elso sorbol szamoljuk ki, ezt pedig hasonlitjuk a tobbi sor r ertekehez
	for (int oszlop_ind = 0; oszlop_ind < this->n; ++oszlop_ind)
	{
		if (this->szomszedsagi_matrix[0][oszlop_ind] != 0)
		{
			r_ertek++;
		}
	}

	bool regularis = true;
	int sor_ind = 1;
	while ((sor_ind < this->n) && regularis)
	{
		int oszlop_ind = 0, jelen_ertek = 0;
		bool lehet_regularis = true;

		while ((oszlop_ind < this->n) && lehet_regularis)
		{
			if (this->szomszedsagi_matrix[sor_ind][oszlop_ind] != 0)
			{
				jelen_ertek++;
				lehet_regularis = (jelen_ertek <= r_ertek);			// Ha mar meghaladtuk az r_ertek-et, akkor felesleges tovabb haladni
			}
			oszlop_ind++;
		}

		regularis = (r_ertek == jelen_ertek);
		sor_ind++;
	}

	return regularis;
}

bool Graf::van_kor() const
{
	vector<bool> latogatott(this->n, false);						// Ha egy csomopont mar latogatott vagy sem
	vector<bool> jelenlegi_ut(this->n, false);						// Egy adott hivasban milyen csomopontok tartoznak a jelenlegi rekurziv fahoz

	bool kor = false;												// Feltetelezzuk, hogy nincs kor a grafban
	int csomopont = 0;
	while ((csomopont < this->n) && !kor)
	{
		if (!latogatott[csomopont])
		{
			kor = this->melysegi_visit_kor(csomopont, latogatott, jelenlegi_ut);
		}
		csomopont++;
	}

	return kor;
}

void Graf::megkeres_utvonal(int ind, int veg, list<El>& utvonal, const vector<int>& parent) const
{
	if (ind != veg)
	{
		this->megkeres_utvonal(parent[ind], veg, utvonal, parent);
		utvonal.push_back(El(parent[ind], ind, this->szomszedsagi_matrix[parent[ind]][ind]));
	}
}

bool Graf::van_negativ_kor(vector<double> tavolsag) const
{
	// Azert, hogy mashol is hasznalhato legyen ez a fuggveny negativ kor felismeresere, ezert ha nem a Bellman-Ford hivja meg, akkor futtatjuk Bellman-Ford algoritmusanak egy reszet
	if (tavolsag.size() == 0)											// Ha meg nem futtattuk Bellman-Ford algoritmusat (egy reszet)
	{
		tavolsag.resize(this->n);
		for (int i = 0; i < this->n; ++i) tavolsag[i] = DBL_MAX;

		// Bellman-Ford algoritmusa
		for (int csomopont = 0; csomopont < this->n - 1; ++csomopont)
		{
			for (const El& el : this->el_lista)
			{
				// Relax minden elre
				if (tavolsag[el.kezd] + el.suly < tavolsag[el.veg])
				{
					tavolsag[el.veg] = tavolsag[el.kezd] + el.suly;
				}
			}
		}
	}

	for (const El& el : el_lista)
	{
		// Relax minden elre
		if (tavolsag[el.kezd] + el.suly < tavolsag[el.veg])				// Ha teljesul, akkor talaltunk egy negativ kort
		{
			return true;
		}
	}
	return false;
}

vector<int> Graf::vegpontok() const
{
	vector<int> vegpontok;

	// Egy csucspont akkor vegpont, ha az incidencia matrixban az adott csucspont soraban csak egyetlen nem 0 ertek van (tehat csak 1 el csatlakozik hozza)
	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		// Csak addig megyunk, mig lehetseges vegpont
		int oszlop_ind = 0, elek_csomopont = 0;

		while ((oszlop_ind < this->m) && elek_csomopont < 2)
		{
			if (this->incidencia_matrix[sor_ind][oszlop_ind] != 0)
			{
				elek_csomopont++;
			}
			oszlop_ind++;
		}

		if (elek_csomopont == 1) vegpontok.push_back(sor_ind);	// Ha vegpont
	}

	return vegpontok;
}

vector<int> Graf::izolalt_csomopontok() const
{
	vector<int> izolalt_pontok;
	for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
	{
		bool izolalt = true;	// eloszor feltetelezzuk, hogy izolalt
		int oszlop_ind = 0;		// atjarjuk a jelenlegi csucs oszlopat

		while ((oszlop_ind < this->n) && izolalt)
		{
			if (this->szomszedsagi_matrix[sor_ind][oszlop_ind] == 0)
			{
				oszlop_ind++;
			}
			else
			{
				izolalt = false;
			}
		}

		if (izolalt) izolalt_pontok.push_back(sor_ind);		// Ha izolalt
	}

	return izolalt_pontok;
}

vector<int> Graf::k_rendu_ismerosok(int ind, int k) const
{
	if (ind < 1 || ind > this->n) throw WrongIndex();
	if (k < 1) throw WrongValue();

	ind--;
	vector<int> tavolsag(this->n, INT32_MAX);
	tavolsag[ind] = 0;

	vector<bool> latogatott(this->n, false);
	latogatott[ind] = true;

	queue<int> sor;
	sor.push(ind);

	int jelen_ut_hossz = 1;									// Eleinte 1 hosszusagu utakat keresunk
	while ((sor.size() != 0) && (jelen_ut_hossz <= k))		// Ha nem lehetseges tovabbi megoldasokat talalni, akkor megallunk
	{
		int jelen_ind = sor.front();
		sor.pop();											// Elso elem torles

		for(const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_ind])
		{
			if (!latogatott[szomszed.ind])
			{
				latogatott[szomszed.ind] = true;
				tavolsag[szomszed.ind] = tavolsag[jelen_ind] + 1;
				jelen_ut_hossz = tavolsag[szomszed.ind];

				sor.push(szomszed.ind);
			}
		}
	}

	vector<int> k_rendu_ismerosok;
	for (int csomopont = 0; csomopont < this->n; ++csomopont)
	{
		if (tavolsag[csomopont] == k) k_rendu_ismerosok.push_back(csomopont);	// Csak ha k tavolsagra van helyezzuk be
	}

	return k_rendu_ismerosok;
}

vector<int> Graf::topologiai_rendezes() const
{
	vector<int> topo;
	vector<bool> latogatott(this->n, false);

	/* 
	 	A topologiai rendezes soran bejarjuk a grafot (Jelen esetben melysegi bejarassal) es egy csomopont akkor kerul a vector - ba
		ha mar feldolgoztuk az osszes szomszedjat. Viszont ez egy forditott topologiai sorrendet ad meg, amit meg kell forditani. 
	*/

	for (int i = 0; i < this->n; ++i)
	{
		if (!latogatott[i])
		{
			this->topologiai_visit(i, latogatott, topo);
		}
	}

	// Megforditjuk a vector-t
	vector<int> topo_vegleges(this->n);
	int q = 0;
	for (int i = static_cast<int>(topo.size() - 1); i >= 0; --i)
	{
		topo_vegleges[q++] = topo[i];
	}

	return topo_vegleges;
}

// Algoritmusok
pair< vector<double>, vector<int> > Graf::Moore_SP(int indulas_ind) const
{
	/*
		Moore (Shortest Path) algoritmusa egy adott sulyozatlan grafban hatarozza meg a legrovidebb tavolsagokat egy csomopontbol
		minden mas tobbi csomopontba.
		Algoritmus mukodese: A indulasi csomopont tavolsaga eleinte 0, ugyanis innen indulva a tavolsagot 0-nak tekinthetjuk.
		Egy varakozasi sort hasznalva behelyezzuk a sorba az indulasi csomopontot es addig haladunk az algoritmussal, amig
		a sor nem ures. Minden egyes iteracioban kiveszunk egy elemet a sorbol es minden szomszedjara megnezzuk, ha meg nem
		latogatott (A tavolsag vegtelen), akkor az adott csomopont tavolsaga a jelenlegi csomopont tavolsaga + 1.
	*/

	indulas_ind--;																	// 0-tol indexelunk
	if (indulas_ind < 0 || indulas_ind >= this->n) throw WrongIndex();				// Helytelen index eseten hibakezeles

	vector<int> parent(this->n, -1);												// Az utak visszakovetesehez szukseges (Egy gyokeres fat alkotnak a minimalis utak, ezert eleg egy vector tarolni)
	vector<double> tavolsag(this->n, DBL_MAX);										// A legrovidebb utat jegyzi meg az 'indulas_ind'-bol (Eleinte mindegyik csomopont INF, tehat elerhetetlen)

	// Varakozasi sor inicializalasa
	queue<int> sor;
	sor.push(indulas_ind);
	tavolsag[indulas_ind] = 0;														// Az indulasi csomopontbol a tavolsag 0

	while (!sor.empty())															// Standard szelessegi bejaras
	{
		int jelen_ind = sor.front();
		sor.pop();

		for (const SzListaElem& szomszed : this->szomszedsagi_lista[jelen_ind])		// Minden szomszed
		{
			if (tavolsag[szomszed.ind] == DBL_MAX)									// Ha meg nem hataroztuk meg a tavolsagot az adott csomopontba
			{
				tavolsag[szomszed.ind] = tavolsag[jelen_ind] + 1;
				parent[szomszed.ind] = jelen_ind;
				sor.push(szomszed.ind);
			}
		}
	}

	return { tavolsag, parent };
}

pair< vector<double>, vector<int> > Graf::Dijkstra_SP(int indulas_ind) const
{
	/*
		Dijkstra (Shortest Path) algoritmusa egy adott sulyozott grafban hatarozza meg a legrovidebb tavolsagokat egy csomopontbol
		minden mas tobbi csomopontba.
		Algoritmus mukodese: A indulasi csomopont tavolsaga eleinte 0, ugyanis innen indulva a tavolsagot 0-nak tekinthetjuk.
		Egy elsobbsegi sort hasznalva mindig azt a csomopontot valasszuk ki, amelynek a legkisebb a tavolsaga az indulasi csomoponttol.
		Bejarjuk a szomszedjait, ha ismeretlenek es tudunk javitani a tavolsagon, akkor behelyezzuk az elsobbseg sorba.
		Az algoritmus ismetlodik, amig az elsobbsegi sor nem ures.
	*/

	indulas_ind--;																	// 0-tol indexelunk
	if (indulas_ind < 0 || indulas_ind >= this->n) throw WrongIndex();				// Helytelen index eseten hibakezeles

	vector<bool> ismeretlen(this->n, true);											// Ismeretlen csomopontok (Hol nem jartunk meg)
	vector<int> parent(this->n, -1);												// Az utak visszakovetesehez szukseges (Egy gyokeres fat alkotnak a minimalis utak, ezert eleg egy vector tarolni)
	vector<double> tavolsag(this->n, DBL_MAX);										// A legrovidebb utat jegyzi meg az 'indulas_ind'-bol (Eleinte mindegyik csomopont INF, tehat elerhetetlen)

	auto Dijkstra_hasonlit = [](const DijkstraNode& elso, const DijkstraNode& masodik) { return elso.tavolsag > masodik.tavolsag; };	// Lambda kifejezes
	priority_queue < DijkstraNode, vector<DijkstraNode>, decltype(Dijkstra_hasonlit) > sor(Dijkstra_hasonlit);

	tavolsag[indulas_ind] = 0;														// Az indulasi csomopontnak a tavolsaga 0
	sor.push( DijkstraNode(indulas_ind, 0) );										// Indulasi index (Ezt fogja valasztani eloszor)

	while (!sor.empty())
	{
		int min_ind = sor.top().csomopont;											// A csomopont amely a legkozelebb van
		double min_tavolsag = sor.top().tavolsag;									// Es a tavolsaga
		sor.pop();

		ismeretlen[min_ind] = false;

		if (tavolsag[min_ind] >= min_tavolsag)										// Optimalizalas: ha mar letezik egy jobb ut az adott csomopontba, akkor nem talalhatunk jobbat (csak akkor lenne lehetseges ha negativ utak is lennenek)
		{
			for (const SzListaElem& szomszed : this->szomszedsagi_lista[min_ind])
			{
				if (ismeretlen[szomszed.ind])										// Ha meg nem latogattuk meg
				{
					double uj_ut = tavolsag[min_ind] + szomszed.suly;				// uj_ut = tavolsag a jelenlegi csomopontba + tavolsag ebbol a csomopontbol a szomszed_csomopont-ba

					if (uj_ut < tavolsag[szomszed.ind])								// Relax ut
					{
						tavolsag[szomszed.ind] = uj_ut;
						sor.push( DijkstraNode(szomszed.ind, uj_ut) );				// Behelyezzuk az elsobbsegi sorba a csomopont indexet es tavolsagat
						// Mivel egy csomopontban tobbszor is frissithetjuk az utat mielott feldolgozzuk, ezert duplikansok is bekerulhetnek (pl: tobbszor kerul be az 1-es csomopont)

						parent[szomszed.ind] = min_ind;								// Megjegyezzuk, hogy hogy jutottunk el a szomszed csomopontba
					}
				}
			}
		}
	}

	return { tavolsag, parent };
}

pair< vector<double>, vector<int> > Graf::Bellman_Ford_SP(int indulas_ind) const
{
	indulas_ind--;
	vector<double> tavolsag(this->n, DBL_MAX);
	tavolsag[indulas_ind] = 0;
	vector<int> parent(this->n, -1);

	for (int csomopont = 0; csomopont < this->n - 1; ++csomopont)
	{
		for (const El& el : this->el_lista)
		{
			// Relax minden elre
			if (tavolsag[el.kezd] + el.suly < tavolsag[el.veg])
			{
				tavolsag[el.veg] = tavolsag[el.kezd] + el.suly;
				parent[el.veg] = el.kezd;
			}
		}
	}

	if (this->van_negativ_kor(tavolsag)) throw GraphHasNegativeCycle();		// Ha van negativ kor a grafban jelezzuk

	return { tavolsag, parent };
}

// << operator
ostream& Graf::kiir(ostream& stream) const
{
	if (this->el_lista.size())
	{
		stream << "El lista: " << '\n';
		for (const El& el : this->el_lista)
		{
			stream << el.kezd + 1 << ' ' << el.veg + 1 << ' ' << el.suly << '\n';
		}
		stream << '\n';
	}

	if (this->szomszedsagi_matrix.size())
	{
		stream << "Szomszedsagi matrix: " << '\n';
		for (int i = 0; i < this->n; ++i)
		{
			stream << i + 1 << ": ";
			for (int j = 0; j < this->n; ++j)
			{
				if (this->szomszedsagi_matrix[i][j] == DBL_MAX)
				{
					stream << "INF ";
				}
				else
				{
					stream << this->szomszedsagi_matrix[i][j] << ' ';
				}
			}
			stream << '\n';
		}
		stream << '\n';
	}

	if (this->incidencia_matrix.size())
	{
		stream << "Incidencia matrix: " << endl;
		for (int i = 0; i < this->n; ++i)
		{
			stream << i + 1 << ": ";
			for (int j = 0; j < this->m; ++j)
			{
				stream << this->incidencia_matrix[i][j] << ' ';
			}
			stream << '\n';
		}
		stream << '\n';
	}

	if (this->szomszedsagi_lista.size())
	{
		stream << "Szomszedsagi lista: " << '\n';
		for (int sor_ind = 0; sor_ind < this->n; ++sor_ind)
		{
			stream << sor_ind + 1 << ": ";

			for (const SzListaElem& szomszed : this->szomszedsagi_lista[sor_ind])
			{
				cout << szomszed.ind + 1 << '(' << szomszed.suly << ") ";
			}
			stream << '\n';
		}
		stream << '\n';
	}

	return stream;
}

ostream& operator<<(ostream& stream, const IranyitatlanGraf& graf)
{
	return graf.kiir(stream);
}

ostream& operator<<(ostream& stream, const IranyitottGraf& graf)
{
	return graf.kiir(stream);
}

// << operator vector
ostream& operator<<(ostream& stream, const list<El>& elek) // cout << elek;
{
	if (elek.size() != 0)
	{
		stream << "Elek: ";
		for(const El& el : elek)
		{
			stream << '(' << el.kezd + 1 << ", " << el.veg + 1 << "){" << el.suly << "} ";
		}
		stream << '\n';
	}

	return stream;
}

ostream& operator<<(ostream& stream, const list<Coord2D>& koordinatak)
{
	stream << "Koordinatak: ";
	for (const Coord2D& koordinata : koordinatak)
	{
		stream << koordinata.x << ',' << koordinata.y << ' ';
	}
	stream << '\n';

	return stream;
}

ostream& operator<<(ostream& stream, const vector<int>& v)
{
	if (v.size() != 0)
	{
		for (const int& elem : v) stream << elem << ' ';
		stream << '\n';
	}
	return stream;
}

ostream& operator<<(ostream& stream, const list<int>& v)
{
	if (v.size() != 0)
	{
		for (const int& elem : v) stream << elem << ' ';
		stream << '\n';
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector<double>& v)
{
	if (v.size() != 0)
	{
		for (int i = 0; i < v.size() - 1; ++i) stream << v[i] << ", ";
		stream << v[v.size() - 1] << '\n';
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector< vector<bool> >& v)
{
	if (v.size() != 0)
	{
		for (int sor = 0; sor < v.size(); ++sor)
		{
			for (int oszlop = 0; oszlop < v[sor].size() - 1; ++oszlop)
			{
				stream << v[sor][oszlop] << ' ';
			}
			if (v[sor].size() != 0)
			{
				stream << v[sor][v[sor].size() - 1] << '\n';
			}
		}
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector< vector<int> >& v)
{
	if (v.size() != 0)
	{
		for (int sor = 0; sor < v.size(); ++sor)
		{
			for (int oszlop = 0; oszlop < v[sor].size() - 1; ++oszlop)
			{
				stream << v[sor][oszlop] << ' ';
			}
			if (v[sor].size() != 0)
			{
				stream << v[sor][v[sor].size() - 1] << '\n';
			}
		}
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector< list<int> >& v)
{
	if (v.size() != 0)
	{
		for (const list<int>& lista : v)
		{
			for (const int& lista_elem : lista)
			{
				stream << lista_elem << ' ';
			}
			stream << '\n';
		}
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector< vector<double> >& v)
{
	if (v.size() != 0)
	{
		for (int sor = 0; sor < v.size(); ++sor)
		{
			for (int oszlop = 0; oszlop < v[sor].size() - 1; ++oszlop)
			{
				stream << v[sor][oszlop] << ' ';
			}
			if (v[sor].size() != 0)
			{
				stream << v[sor][v[sor].size() - 1] << '\n';
			}
		}
	}
	return stream;
}

ostream& operator<<(ostream& stream, const vector< vector<Coord2D> >& v)
{
	if (v.size() != 0)
	{
		for (int sor = 0; sor < v.size(); ++sor)
		{
			for (int oszlop = 0; oszlop < v[sor].size(); ++oszlop)
			{
				stream << v[sor][oszlop].x << ',' << v[sor][oszlop].y << ' ';
			}
			stream << '\n';
		}
	}
	return stream;
}
/*-------------------------------------------------------------------- Graf ---------------------------------------------------------------------------------*/
