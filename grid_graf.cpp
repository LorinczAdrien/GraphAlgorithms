#include "graf.h"

/*------------------------------------------------------------------- Grid Graf -----------------------------------------------------------------------------*/
// Coord2D
ostream& operator<<(ostream& stream, const vector<Coord2D>& koordinatak)
{
	for (const Coord2D& koordinata : koordinatak)
	{
		stream << "x = " << koordinata.x << ", y = " << koordinata.y << '\n';
	}

	return stream;
}

// Class
GridGraph::GridGraph(string file_in)
{
	ifstream fin(file_in);

	if (fin.fail()) throw FileNotAccesible();

	fin >> this->sor >> this->oszlop;
	this->grid.resize(this->sor);
	for (int i = 0; i < this->sor; ++i) this->grid[i].resize(this->oszlop);

	// Beolvas ertekek
	for (int i = 0; i < this->sor; ++i)
	{
		for (int j = 0; j < this->oszlop; ++j)
		{
			fin >> this->grid[i][j];
		}
	}
}

GridGraph::GridGraph(int sor, int oszlop, const vector< vector<int> >& adott_grid)
{
	this->sor = sor; this->oszlop = oszlop;
	this->grid = adott_grid;
}

// Operators
GridGraph& GridGraph::operator =(const GridGraph& grid_graf)
{
	if (this != &grid_graf)
	{
		this->sor = grid_graf.sor;
		this->oszlop = grid_graf.oszlop;
		this->grid = grid_graf.grid;
	}
	return *this;
}

// Utility
list<Coord2D> GridGraph::get_szomszedok(const Coord2D& jelen_poz, int mode) const
{
	if (sor < 0 || oszlop < 0 || jelen_poz.x >= this->sor || jelen_poz.y >= this->oszlop) throw IndexOutOfBounds();

	list<Coord2D> szomszedok;											// Visszateritendo lista
	vector<Coord2D> szomszed_lepesek{ Coord2D(-1, 0), Coord2D(0, 1), Coord2D(1, 0), Coord2D(0, -1), Coord2D(-1, 1), Coord2D(1, 1), Coord2D(1, -1), Coord2D(-1, -1) };

	if (this->grid[jelen_poz.x][jelen_poz.y] == 0) return szomszedok;	// Ha fal van ott akkor nem teritunk szomszedokat

	for (int szomszed = 0; szomszed < mode; ++szomszed)
	{
		if (is_valid_coordinate(jelen_poz + szomszed_lepesek[szomszed]))
		{
			szomszedok.push_back(jelen_poz + szomszed_lepesek[szomszed]);
		}
	}
	return szomszedok;
}

list<Coord2D> GridGraph::get_szabad_szomszedok(const Coord2D& jelen_poz, int mode) const
{
	if (sor < 0 || oszlop < 0 || jelen_poz.x >= this->sor || jelen_poz.y >= this->oszlop) throw IndexOutOfBounds();

	list<Coord2D> szomszedok;											// Visszateritendo lista
	vector<Coord2D> szomszed_lepesek{ Coord2D(-1, 0), Coord2D(0, 1), Coord2D(1, 0), Coord2D(0, -1), Coord2D(-1, 1), Coord2D(1, 1), Coord2D(1, -1), Coord2D(-1, -1) };

	if (this->grid[jelen_poz.x][jelen_poz.y] == 0) return szomszedok;	// Ha fal van ott akkor nem teritunk szomszedokat

	for (int szomszed = 0; szomszed < mode; ++szomszed)
	{
		Coord2D szomszed_poz = jelen_poz + szomszed_lepesek[szomszed];
		if (is_valid_coordinate(szomszed_poz) && (this->grid[szomszed_poz.x][szomszed_poz.y] == 1))
		{
			szomszedok.push_back(szomszed_poz);
		}
	}
	return szomszedok;
}

void GridGraph::megkeres_ut(Coord2D jelen_poz, list<Coord2D>& ut, const Coord2D& keresett, const vector< vector<Coord2D> >& parent) const
{
	if (jelen_poz != keresett)
	{
		this->megkeres_ut(parent[jelen_poz.x][jelen_poz.y], ut, keresett, parent);
		ut.push_back(jelen_poz);
	}
}

bool GridGraph::is_valid_coordinate(const Coord2D& kord) const
{
	return ((kord.x >= 0) && (kord.y >= 0) && (kord.x < this->sor) && (kord.y < this->oszlop));
}

Coord2D GridGraph::D1_to_D2(int poz) const
{
	Coord2D D2_koordinata(0, 0);
	
	while (poz >= this->oszlop)		// Hanyadik sor
	{
		D2_koordinata.x++;
		poz -= this->oszlop;
	}
	D2_koordinata.y = poz;			// A maradek az oszlop index

	return D2_koordinata;
}

int GridGraph::D2_to_D1(const Coord2D& poz) const
{
	return poz.x * this->oszlop + poz.y;
}

// << operator
ostream& GridGraph::kiir(ostream& stream) const
{
	stream << "Grid: " << this->sor << 'x' << this->oszlop << '\n';
	for (int i = 0; i < this->grid.size(); ++i)
	{
		for (int j = 0; j < this->grid[i].size(); ++j)
		{
			stream << ' ' << this->grid[i][j] << ' ';
		}
		stream << '\n';
	}
	stream << '\n';
	return stream;
}

ostream& operator <<(ostream& stream, const GridGraph& graf)
{
	return graf.kiir(stream);
}

// User functions
void GridGraph::set_ertek(int sor, int oszlop, int uj_ertek)
{
	if (sor < 0 || oszlop < 0 || sor >= this->sor || oszlop >= this->oszlop) throw IndexOutOfBounds();
	this->grid[sor][oszlop] = uj_ertek;
}

int GridGraph::get_ertek(int sor, int oszlop) const
{
	if (sor < 0 || oszlop < 0 || sor >= this->sor || oszlop >= this->oszlop) throw IndexOutOfBounds();
	return this->grid[sor][oszlop];
}

list<int> GridGraph::get_szabad_szomszedok_1D(int poz, int mode) const
{
	list<int> szomszedok_1D;
	Coord2D jelen_poz = this->D1_to_D2(poz);

	list<Coord2D> szomszedok_2D = this->get_szabad_szomszedok(jelen_poz, 4);
	for (const Coord2D& szomszed : szomszedok_2D)
	{
		szomszedok_1D.push_back(this->D2_to_D1(szomszed));
	}

	return szomszedok_1D;
}

pair< list<Coord2D>, int > GridGraph::BFS_SP(const Coord2D& start, const Coord2D& end) const
{
	if (!is_valid_coordinate(start) || !is_valid_coordinate(end)) throw IndexOutOfBounds();					// Ha nem 'valid' koordinatak, vagyis ha kint vannak a 'grid'-en
	if ((this->grid[start.x][start.y] == 0) || (this->grid[end.x][end.y] == 0)) throw WallAtCoordinate();	// Ha fal van az adott koordinataknal

	list<Coord2D> ut;															// Visszateritett ut
	int ut_hossz = 0;															// Es hossza

	bool megtalalt = false;
	int jelen_reteg_maradt = 1, kovetkezo_reteg = 0;
	queue<Coord2D> sor;
	vector< vector<bool> > latogatott(this->sor, vector<bool>(this->oszlop, false));
	vector< vector<Coord2D> > parent(this->sor, vector<Coord2D>(this->oszlop, Coord2D(-1, -1)));

	// Eloszor is behelyezzuk az indulasi koordinatakat
	sor.push(start);
	latogatott[start.x][start.y] = true;										// Mostmar latogatott

	while (!sor.empty())														// Ameddig meg tudunk keresni es nem talaltuk meg, addig keres
	{
		Coord2D jelen_poz = sor.front();										// Jelenlegi pozicio
		sor.pop();

		if (jelen_poz == end)
		{
			megtalalt = true;
			break;
		}

		list<Coord2D> szomszedok = this->get_szabad_szomszedok(jelen_poz, 8);
		for (const Coord2D& szomszed_poz : szomszedok)
		{
			if (!latogatott[szomszed_poz.x][szomszed_poz.y])					// Ha meg nem jartunk azon a szabad helyen
			{
				latogatott[szomszed_poz.x][szomszed_poz.y] = true;
				sor.push(szomszed_poz);
				kovetkezo_reteg++;

				parent[szomszed_poz.x][szomszed_poz.y] = jelen_poz;				// Megjegyezzuk, hogy honnan jottunk
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
		ut.push_back(start);
		this->megkeres_ut(end, ut, start, parent);

		vector< vector<bool> > uton(this->sor, vector<bool>(this->oszlop, false));
		for (const Coord2D& kord : ut) uton[kord.x][kord.y] = true;

		cout << "Megoldas:\n";
		for (int i = 0; i < this->sor; ++i)
		{
			for (int j = 0; j < this->oszlop; ++j)
			{
				if (uton[i][j]) cout << " x ";
				else cout << " 0 ";
			}
			cout << '\n';
		}
		cout << '\n';
	}
	else
	{
		cout << "Nincs megoldas: Nem elerheto az end point!";
		ut_hossz = -1;
	}

	return { ut, ut_hossz };
}

pair< list<Coord2D>, int > GridGraph::A_star_SP(const Coord2D& start, const Coord2D& end) const
{
	if (!is_valid_coordinate(start) || !is_valid_coordinate(end)) throw IndexOutOfBounds();					// Ha nem 'valid' koordinatak, vagyis ha kint vannak a 'grid'-en
	if ((this->grid[start.x][start.y] == 0) || (this->grid[end.x][end.y] == 0)) throw WallAtCoordinate();	// Ha fal van az adott koordinataknal

	list<Coord2D> ut;															// Visszateritett ut
	int ut_hossz = 0;															// Es hossza
	bool vegpont_elerheto = false;

	list<AStarNode> open, closed;
	open.push_back( AStarNode(start, start, start, end) );						// A nyitott szetbe behelyezzuk az indulasi csomopontot

	while (!open.empty())
	{
		AStarNode jelen_node;
		jelen_node.total_distance = DBL_MAX;

		auto it = open.begin(); auto legkisebb = it;
		while (it != open.end())												// Megkeressuk a minimalis ossz tavolsagu csomopontot
		{
			if (it->total_distance < jelen_node.total_distance)
			{
				legkisebb = it;													// Lementjuk a poziciot, a kesobbi torles erdekeben
				jelen_node.total_distance = it->total_distance;					// Frissit minimum												
			}
			it++;
		}

		// Kivesszuk a nyitott listabol es behelyezzuk a zart-ba
		jelen_node = *legkisebb;
		open.erase(legkisebb);
		closed.push_back(jelen_node);

		// Ha elertunk a celhoz akkor kilephetunk
		if (jelen_node.poz == end)
		{
			vegpont_elerheto = true;
			break;
		}

		list<Coord2D> szomszedok = this->get_szabad_szomszedok(jelen_node.poz, 8);
		for (const Coord2D& szomszed : szomszedok)
		{
			// Ha a zart listaban van akkor lephetunk tovabb
			bool szomszed_zart = false;
			for (const AStarNode& zart_node : closed)
			{
				if (szomszed == zart_node.poz)
				{
					szomszed_zart = true;
					break;
				}
			}
			if (szomszed_zart) continue;

			// Nincs a zart listaban
			// Megnezzuk, ha meg nincs a nyilt listaban akkor behelyezzuk, ellenkezo esetben pedig, ha tudunk javitani az ut koltsegen akkor javitunk
			bool szomszed_nyilt = false;
			auto szomszed_node = open.begin();
			for (auto nyilt_node = open.begin(); nyilt_node != open.end(); nyilt_node++)
			{
				if (szomszed == nyilt_node->poz)
				{
					szomszed_nyilt = true;
					szomszed_node = nyilt_node;
					break;
				}
			}

			if (szomszed_nyilt)													// Mar a listaban, ezert megnezzuk, ha tudunk javitani rajta
			{
				double uj_tavolsag = jelen_node.local_distance + jelen_node.poz.tavolsag_ket_koordinata_kozott(szomszed_node->poz) + szomszed_node->poz.tavolsag_ket_koordinata_kozott(end);
				if (uj_tavolsag < szomszed_node->total_distance)
				{
					szomszed_node->total_distance = uj_tavolsag;
					szomszed_node->local_distance = jelen_node.local_distance + jelen_node.poz.tavolsag_ket_koordinata_kozott(szomszed_node->poz);
				}
			}
			else																// Nincs a listaban, ezert behelyezzuk
			{
				open.push_back( AStarNode(szomszed, jelen_node.poz, start, end) );
			}
		}
	}

	if (vegpont_elerheto)
	{
		// Visszakeressuk az utvonalat
		Coord2D jelen_node(end); ut.push_back(end);
		while (jelen_node != start)
		{
			for (const AStarNode& node : closed)
			{
				if (node.poz == jelen_node)
				{
					jelen_node = node.parent_poz;
					ut.push_back(jelen_node);
					ut_hossz++;
				}
			}
		}
		ut.reverse();															// Megforditjuk a listat, ugyanis a vegetol kerestuk az utat, igy egy forditott utvonalat kaptunk
	}

	return { ut, ut_hossz };
}

/*------------------------------------------------------------------- Grid Graf -------------------------------------------------------------------------------------*/