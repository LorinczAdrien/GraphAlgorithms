#include "graf.h"

int main()
{
	GridGraph graf("grid_in.txt");
	list<Coord2D> ut; int ut_hossz = 0;

	Coord2D start(0, 0), end(0, 8);
	tie(ut, ut_hossz) = graf.BFS_SP(start, end);
	cout << "Ut: " << ut;
	cout << "Ut hossza: " << ut_hossz << '\n';

	tie(ut, ut_hossz) = graf.A_star_SP(start, end);
	cout << "Ut: " << ut;
	cout << "Ut hossza: " << ut_hossz << '\n';

	return 0;
}