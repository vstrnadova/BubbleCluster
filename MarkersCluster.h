/* Copyright 2014 Veronika Strnadova-Neeley, Aydin Buluc, Jarrod Chapman
 * 
 *  MarkersCluster.h is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MarkersCluster.h is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MarkersCluster.h.  If not, see <http://www.gnu.org/licenses/>.  
 */
#ifndef _MARKERS_CLUSTER_H
#define _MARKERS_CLUSTER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <functional>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <list>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <sys/time.h>
#include "MCHelpers.h"

using namespace std;


class MarkersCluster
{
public:
    MarkersCluster() {};
    MarkersCluster(ifstream & input1, char* argv[], double sens);
    tuple<double,double,double> returnlod(int i, int j);
    int mergerepptlists(vector<list<int>> & repptsvector, int i,int j, double threshold, double maxlod);
    void assignlqmarkers(list<int> & lqmarkers, vector<list<int> > & clusters, vector<list<int> > & reppts, int * markerassignedarray, double lowqualc);
    void checkBoundary(int clusternum, std::list<int> & reppts, int highestreppt, int marker, double maxlod, double threshold);
    void bubbleCluster(vector<list<int>> & clustervec, vector<list<int>> & repptvec, int *randommarkerarray,int nmarkers,
                       int *markerassigned, double threshold, double selflodthreshold);
    void mergesmallclusters(vector<list<int>> & clusters, vector<list<int>> & repptlists, double sigma, double lowqualc, double threshold);
    void printclusters(vector<list<int>> & clusters, char* ofname);
    void printreppoints(vector<list<int>> & reppointvec, char* ofname);
    vector<int> getreppointvec(vector<list<int>> & reppointvec, int index);
    int NMarkers() { return nmarkers; };
    int Nindividuals() { return nindividuals; };
    int getNumNonmissing(int markernum);
    void printMarkerRptAssignments(int *markerassignedarray, int nummarkers);
    void printRptLODMtx(list<int> & rpptlist, int cnum); 
private:
    vector< vector<sparsecall> > calls1;
    vector< vector<sparsecall> > calls2;
    vector<string> markers1; // names of markers are in markers1 and markers2
    vector<string> markers2;
    int nmarkers;
    int nindividuals;
    double sensitivity;
    double c;
    double sigma;
};
    
#include "MarkersCluster.cpp"
#endif
