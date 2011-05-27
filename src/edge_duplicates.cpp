/* Function for finding edge loops, duplicates and bi-directional edges.
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>

extern "C" {
#include <R.h>

using namespace std;


void edgeDuplicates(int *edgeA, int *edgeB, int *numedg, int *loops, int *dups, bool *verbose) 

	{

	int i, j;
	float prog;
	bool add = TRUE;
	set<int> tset;
	set<int> nset;
	set<int> common;
	vector< set<int> > esets; // Edge sets.

	// Loop through edges and ask if they're loops, duplicates, or bidirectional as we go by assembling a vector of sets.
	for(i = 0; i < *numedg; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numedg-1)*100;

			Rprintf("\r   Checking for loops and duplicate edges... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(edgeA[i] == edgeB[i]){
			loops[i] = 1;
			continue;
		}else{
			tset.insert(edgeA[i]);
			tset.insert(edgeB[i]);

			for(j = 0; j < esets.size(); j++){

				nset = esets.at(j);
				set_intersection(tset.begin(),tset.end(), nset.begin(),nset.end(), inserter(common, common.begin()));

				if(common.size() == 2){ // Duplicate or bi-directional edge.
					dups[i] = 1;
					nset.clear();
					common.clear();
					add = FALSE;
					break;
					}

				nset.clear();
				common.clear();
				}

			if(add){
				esets.push_back(tset);
				}
			tset.clear();
			add = TRUE;

			}

		}

	}


	}





