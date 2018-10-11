/* Copyright 2014 Veronika Strnadova-Neeley, Aydin Buluc, Jarrod Chapman
 * 
 *  BubbleClusterV4.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  BubbleClusterV4.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BubbleClusterV4.cpp.  If not, see <http://www.gnu.org/licenses/>.  
 */
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
#include "MarkersCluster.h"

using namespace std;
//set MERGESMALLCLUSTERS to 0 if you do NOT want to attempt to merge small clusters
// (defined by clusters of size < smallclustersize) with larger clusters
const int MERGESMALLCLUSTERS=1;
MarkersCluster MC;

bool selfLOD (int i,int j) {
    return (get<2>(MC.returnlod(i,i)) > get<2>(MC.returnlod(j,j)));
}

bool numNonmissing(int i, int j){
   return (MC.getNumNonmissing(i) > MC.getNumNonmissing(j));
}

int main(int argc, char* argv[])
{
    if(argc != 7)
    {
		cout << "Usage: ./v4 <markerinputfile> <lodthreshold> <nonmissingthreshold> <c> <smallclustersize> <output>" << endl;
		cout << "The input file should have one marker name per line, with tab-separated"
		     << " genotypes following each marker name."<<endl;
		return 0;
    }
    ifstream input1(argv[1], ifstream::in);

    // variables for timing purposes
    timeval TclusterS, TclusterF;
    timeval TnmsortS, TnmsortF;
    timeval TlowqualS, TlowqualF;
    timeval TsmallcsS, TsmallcsF;
	
    MC = MarkersCluster(input1, argv, 0.5);
    int nmarkers = MC.NMarkers();
	
    auto start = std::chrono::system_clock::now();
		
    //initialize an array of markers, which keeps track of which marker has been assigned to a group
    // markerassigned[] is going to hold the rep pt that each marker is assigned to
    // int markerassigned[nmarkers];
    int markersarray[nmarkers];
    int markerassigned[nmarkers];
    for(int i = 0; i<nmarkers; i++)
    {
      markersarray[i] = i;
    }
 	
    cout << "created markersarray, markerassigned"<<endl;	
    vector<int> markersvector (markersarray, markersarray+nmarkers);

    // sort markers by number of non-missing entries
    gettimeofday(&TnmsortS,NULL);
       sort(markersvector.begin(), markersvector.end(), numNonmissing);
    gettimeofday(&TnmsortF,NULL);
	
    int medianpos = markersvector.size()/2;
    int medianmarker = markersvector.at(medianpos);
    cout << "median number of non-missing entries = " << MC.getNumNonmissing(medianmarker)  << endl;
    cout << "median self-LOD = " <<get<2>(MC.returnlod(medianmarker,medianmarker))  << endl;
    long seconds  = TnmsortF.tv_sec  - TnmsortS.tv_sec;
    long useconds = TnmsortF.tv_usec - TnmsortS.tv_usec;
    double nmsort_time = seconds + useconds/1000000.0;
	
    double threshold = atof(argv[2]);
    cout << "threshold = "<<threshold<<endl;
    int nonmissingthreshold = atoi(argv[3]);
    double selflodthreshold = log10(2)*nonmissingthreshold;	       
    cout << "non-missing threshold = "<<nonmissingthreshold<<endl;
    cout << "self-lod threshold = "<< selflodthreshold <<endl;
    double lowqualc = atof(argv[4]); 
    cout << "c = "<<lowqualc<<endl;
    double sigma = atof(argv[5]);
    cout << "sigma = "<<sigma<<endl;
	 
    for(int i=0; i<nmarkers; i++){
       markersarray[i] = markersvector.at(i);
    }
        
    for(int i=0; i<nmarkers;i++){
       if(get<2>(MC.returnlod(markersarray[i],markersarray[i])) < selflodthreshold){
	  cout << "self-LOD less than selflodthreshold for "<<(nmarkers-1-i) << " markers"<<endl; 
	  break;
       }
    }

    std::list<std::list<int>> clusters;
    std::list<std::list<int>> listofrppts;
    list<int> clusterarr[15000];
    list<int> rpptarr[15000];

    cout << "initialized arrays of cluster lists and sketch point lists"<<endl;	
    vector<list<int>> clustvec (clusterarr, clusterarr+sizeof(clusterarr)/sizeof(list<int>));
    vector<list<int>> repvec (rpptarr, rpptarr+sizeof(rpptarr)/sizeof(list<int>));
    // use bubbleCluster to cluster the markers 
    cout << "starting bubbleCluster()..."<<endl;
    gettimeofday (&TclusterS, NULL);
      MC.bubbleCluster(clustvec, repvec, markersarray,nmarkers,markerassigned,threshold,selflodthreshold);
    gettimeofday (&TclusterF, NULL);
    seconds  = TclusterF.tv_sec  - TclusterS.tv_sec;
    useconds = TclusterF.tv_usec - TclusterS.tv_usec;
    double cluster_time = seconds + useconds/1000000.0;
	
    //create a list of low-quality markers
    cout << "initial clustering stage finished; now adding low-quality markers" << endl;
    int mcount = 0;
    list<int> lowqualmarkers;
    while(mcount<nmarkers)
    {
       	   if(get<2>(MC.returnlod(markersarray[mcount],markersarray[mcount])) <= selflodthreshold)
	   {
	      lowqualmarkers.push_back(markersarray[mcount]);	
	   }
	   mcount++;
    }	
    cout << "lowqualmarkers.size() before assignlqmarkers = "<<lowqualmarkers.size()<<endl; 
	
    // assign the low-quality markers to their best-fit cluster
    gettimeofday(&TlowqualS,NULL);
       MC.assignlqmarkers(lowqualmarkers, clustvec, repvec, markerassigned, lowqualc);
    gettimeofday(&TlowqualF,NULL);
    cout << "lowqualmarkers.size() after assignlqmarkers = "<<lowqualmarkers.size()<<endl; 
    lowqualmarkers.clear();
    seconds  = TlowqualF.tv_sec  - TlowqualS.tv_sec;
    useconds = TlowqualF.tv_usec - TlowqualS.tv_usec;
    double lowqual_time = seconds + useconds/1000000.0;

    // report cluster sizes if greater than 1 
    int total=0;	
    for(int ccount=0; ccount<clustvec.size(); ccount++){
	if((clustvec[ccount]).size() > 1)
	      cout << "cluster " << ccount<< " size = " << (clustvec[ccount]).size()<<endl;
        total += (clustvec[ccount]).size();
    }
    cout << "total number of markers in clusters = "<<total<<endl;        

    // attempt to merge small clusters with large ones
    double smallcs_time;
    if(MERGESMALLCLUSTERS){
       gettimeofday(&TsmallcsS,NULL);
          MC.mergesmallclusters(clustvec,repvec,sigma,lowqualc,threshold);
       gettimeofday(&TsmallcsF,NULL);
       seconds  = TsmallcsF.tv_sec  - TsmallcsS.tv_sec;
       useconds = TsmallcsF.tv_usec - TsmallcsS.tv_usec;
       smallcs_time = seconds + useconds/1000000.0;
    }
    
    //print out each cluster and each rppt list to a file
    MC.printclusters(clustvec, argv[6]); 
    MC.printreppoints(repvec, argv[6]);
    int nummarkers = MC.NMarkers();
    
    printf("threshold LOD = %f\n",threshold);
    cout << "non-missing entry sort finished in: " << nmsort_time << " seconds. "<< endl;
    cout << "clustering finished in: " << cluster_time << " seconds. "<< endl;
    cout << "low quality marker placement finished in: " << lowqual_time << " seconds. "<< endl;
    if(MERGESMALLCLUSTERS) {
	cout << "post-processing small clusters finished in: " << smallcs_time << " seconds. "<< endl;
    } else {
        cout << "no small cluster merging "<<endl;
    }
    cout << "total time: non-missing sort+cluster+lq placement+smallclust = " << nmsort_time+cluster_time+lowqual_time+smallcs_time<<endl;
    return 0;
} 
