/* Function for calculating link densities and the partition density for clusters generated within "getLinkCommunities".
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R.
 * Based on the Link communities derived from the algorithm in:
 * Ahn et al. (2010). Link communities reveal multiscale complexity in networks. Nature 466:761-765.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <set>

extern "C" {
#include <R.h>

using namespace std;

void getLinkDensities(int *ma, int *mb, int *ea, int *eb, int *numedg, int *clusnums, double *pdens, double *heights, double *pdmax, int *csize, bool *removetrivial, 
			bool *carriageret, bool *verbose)

	{

	int i, j, k, p = 0, ne, nn, count = 0, csum = clusnums[0], perc = 0;
	float prog;
	double ldens, maxp, best = 0;
	vector<int> mergeA;
	vector<int> mergeB;
	vector<int> edgeA;
	vector<int> edgeB;
	vector<int> drop;
	vector< vector<int> > clusters;
	vector<int> current;
	vector<int> currentTemp;
	vector<int> bestC;
	set<int> nodes;

	copy(ma, ma + *numedg-1, back_inserter(mergeA));
	copy(mb, mb + *numedg-1, back_inserter(mergeB));

	copy(ea, ea + *numedg, back_inserter(edgeA));
	copy(eb, eb + *numedg, back_inserter(edgeB));

	remove("linkcomm_clusters.txt");

	ofstream outfile;
	outfile.open("linkcomm_clusters.txt", ios::out);
	if(! outfile.is_open()){
			Rprintf("\nERROR: can't open linkcomm_clusters.txt!\n"); return;
			}

	for(i = 0; i < (*numedg-1); i++){ // Construct matrix for clusters.
			clusters.push_back(vector<int>());
			}

	// Loop through merges from lowest to highest.
	// At each merge extract cluster elements and write to disk.
	for(i = 0; i < *numedg-1; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numedg-2)*100;

			if(*carriageret){
				Rprintf("   Calculating link densities... %3.2f%\r",prog);
			}else{
				if(i == 0){
					Rprintf("   Calculating link densities...\r\n");
					}
				if(prog >= perc){
					Rprintf("|");
					perc++;
					}
				}

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(mergeA.at(i) < 0 && mergeB.at(i) < 0){
			clusters[i].push_back(-1*mergeA.at(i));
			clusters[i].push_back(-1*mergeB.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i); // Update current cluster subs.
		}else if(mergeA.at(i) > 0 && mergeB.at(i) < 0){
			clusters[i] = clusters[mergeA.at(i)-1];
			clusters[i].push_back(-1*mergeB.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
		}else if(mergeA.at(i) < 0 && mergeB.at(i) > 0){
			clusters[i] = clusters[mergeB.at(i)-1];
			clusters[i].push_back(-1*mergeA.at(i));
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeB.at(i) - 1);
		}else{
			clusters[i].insert(clusters[i].end(), clusters[mergeA.at(i)-1].begin(), clusters[mergeA.at(i)-1].end());
			clusters[i].insert(clusters[i].end(), clusters[mergeB.at(i)-1].begin(), clusters[mergeB.at(i)-1].end());
			sort(clusters[i].begin(),clusters[i].end());
			current.push_back(i);
			drop.push_back(mergeA.at(i) - 1);
			drop.push_back(mergeB.at(i) - 1);
			}

		if(i+1 == csum){ // Reached end of this height, calculate link densities for the current clusters.
			currentTemp = current;
			sort(drop.begin(),drop.end());
			for(j = 0; j < drop.size(); j++){
				currentTemp.erase(currentTemp.begin() + drop.at(j) - j);
				}
			ldens = 0;
			for(j = 0; j < currentTemp.size(); j++){
				// Number of edges.
				ne = clusters[currentTemp.at(j)].size();
				// Number of nodes.
				for(k = 0; k < ne; k++){
					nodes.insert(edgeA.at(clusters[currentTemp.at(j)].at(k)-1));
					nodes.insert(edgeB.at(clusters[currentTemp.at(j)].at(k)-1));
					}
				nn = nodes.size();
				ldens = ldens + (double(ne)*(double(ne)-double(nn)+1.0))/((double(nn)-2.0)*(double(nn)-1.0));
				nodes.clear();
				}
			
			pdens[p] = (2.0/ *numedg)*ldens;
			maxp = pdens[p];

			if(maxp > best){
				best = maxp;
				*pdmax = heights[p];
				bestC = currentTemp;
				}

			p++;
			csum = csum + clusnums[p];
			}

		}

	// Write optimal clusters to disk
	for(i = 0; i < bestC.size(); i++){
		if(*removetrivial && clusters[bestC.at(i)].size() == 2){count++; continue;} // Delete trivial clusters of size 2.
		for(j = 0; j < clusters[bestC.at(i)].size(); j++){
			outfile << clusters[bestC.at(i)].at(j) << " ";
			}
		outfile << endl;
		}

	*csize = bestC.size() - count;
	
	outfile.close();

	}

	}






	
