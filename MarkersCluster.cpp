/* Copyright 2014 Veronika Strnadova-Neeley, Aydin Buluc, Jarrod Chapman
 * 
 *  MarkersCluster.cpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MarkersCluster.cpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MarkersCluster.cpp.  If not, see <http://www.gnu.org/licenses/>.  
 */
tuple<double,double,double> MarkersCluster::returnlod(int i, int j)
{
        
	pair<int,int> flops_accum = intersect_accumulate(calls1[i].begin(), calls1[i].end(), calls2[j].begin(), calls2[j].end(), multcalls());
	// nSame + nOpp = flops
	// nSame - nOpp = accum
	// nSame = (flops + accum) / 2
	int flops = flops_accum.first;
	int accum = flops_accum.second;
	double nSame = (flops+accum) / 2;
	double nOpp = nSame - accum;
	
	double p = nSame / flops; // flops is nInformative
	
	if (p <= sensitivity || p > (1-sensitivity))
	{
		double r = std::min(nSame, nOpp);
		double nr = std::max(nSame, nOpp);
		
		//if(r==0) cout << "r==0"<<endl;
		// [ABAB: is it the case that recombination rate can't be zero]
		//r = (r == 0) ? 0.5: r; // Is this the right choice of pseudo-counts?
		//if(r==0) cout << "after check: r==0"<<endl;
		double theta = r/ (r+nr);
		//if(theta==0) return make_tuple(nSame,nOpp,std::log10(pow(2,r)));
                double lod;
                if(r!=0){
  	         lod = std::log10 ( ( ( pow((1-theta),nr)  * pow(theta,r) ) / pow(0.5,r+nr) ) );
		} else {
  		 lod = nr*std::log10(2);
		}
	 	return make_tuple(nSame,nOpp,lod);
	}
	return make_tuple(nSame,nOpp, (double) 0.0);
}

bool greaterclustersize(pair<int,int> a, pair<int,int> b){
  return (a.second > b.second);
}

bool greaterthan (int i,int j) { return (i>j); }

/* mergerepptlists(repptsvector.at(i), repptsvector.at(j)):
 *   merges two lists of representative points, according to the lod
 *   scores between the boundaries of the points. If unable to merge
 *   based on LOD scores, returns an empty list
 * */
int MarkersCluster::mergerepptlists(vector<list<int>> & repptsvector, int i,int j, double threshold, double maxlod)
{
    int ok = 0;
    
    // here we assume that m is between the boundaries of reppt lists i and j 
    int firstreppta = repptsvector[i].front();
    int lastreppta = repptsvector[i].back();
    int firstrepptb = repptsvector[j].front();
    int lastrepptb = repptsvector[j].back();
    
    double a1b1 = get<2>(returnlod(firstreppta,firstrepptb));
    double a1b2 = get<2>(returnlod(firstreppta,lastrepptb));
    double a2b1 = get<2>(returnlod(lastreppta,firstrepptb));
    double a2b2 = get<2>(returnlod(lastreppta,lastrepptb));
    
    if(repptsvector[i].size()==0 || repptsvector[j].size()==0){
        return -1;
    }
    
    if(repptsvector[i].size()==1 || repptsvector[j].size()==1){
        //cout << " one of the lists' size is 1"<<endl;
        if(repptsvector[i].size() == 1 && repptsvector[j].size() == 1){
            //cout << "in mergerepptlists: both lists are only one point"<<endl;
            repptsvector[i].splice(repptsvector[i].end(),repptsvector[j]);
            return 0;
        } else if(repptsvector[i].size()==1){
            //cout << "in mergerepptlists: reppt list "<<i<<" is only one point"<<endl;
            if(a1b1>=a1b2){
                //cout << "adding reppoint to front of list "<<j<< endl;
                repptsvector[j].push_front(repptsvector[i].front());
                repptsvector[i].clear();
                return 0;
            } else if(a1b2>= a1b1) {
                //cout << "adding reppoint to back of list "<<j<<endl;
                repptsvector[j].push_back(repptsvector[i].front());
                repptsvector[i].clear();
                return 0;
            } else {
                //cout << " could not determine where to place single reppoint "<< repptsvector[i].front()<< " in  repptlist "<<i<<endl;
                repptsvector[i].splice(repptsvector[i].end(),repptsvector[j]);
                return 0;
            }
        } else if (repptsvector[j].size() == 1){
            //cout << "in mergerepptlists: reppt list "<<j<<" is only one point"<<endl;
            if(a1b1>=a2b1){
                //cout << "adding reppoint to front of list "<<i<<endl;
                repptsvector[i].push_front(repptsvector[j].front());
                repptsvector[j].clear();
                return 0;
            } else if (a2b1>=a1b1) {
                //cout << "adding reppoint to back of list "<<i<<endl;
                repptsvector[i].push_back(repptsvector[j].front());
                repptsvector[j].clear();
                return 0;
            } else {
                //cout << " could not determine where to place single reppoint "<< repptsvector[j].front()<< " in  repptlist "<<j<<endl;
                repptsvector[i].splice(repptsvector[i].end(),repptsvector[j]);
                return 0;
                //return -1;
            }
        }
    }
    // if both rppt lists are of size >1, merge them based on the LOD
    // score between their boundary points a1, a2 and b1, b2	
    if((a1b1 <= a1b2) && (a2b2 <= a1b2)){
            repptsvector[i].reverse();
            repptsvector[j].reverse();
      	    repptsvector[i].splice(repptsvector[i].end(),repptsvector[j]);
            return ok;
        } else if((a1b1 >= a1b2) && (a1b1 >= a2b1) ) {
            repptsvector[i].reverse();
            repptsvector[i].splice(repptsvector[i].end(), repptsvector[j]);
            return ok;
        } else if ((a1b2 <= a2b2) && (a2b1 <= a2b2)){
            repptsvector[j].reverse();
            repptsvector[i].splice(repptsvector[i].end(), repptsvector[j]);
            return ok;
        } else if ((a1b1 <= a2b1) && (a2b2 <= a2b1)){
            repptsvector[i].splice(repptsvector[i].end(), repptsvector[j]);     	    
            return ok;
        } else {
	    cout << "could not determine correct order to merge clusters "<<i<<" and "<<j<<endl;
            repptsvector[i].splice(repptsvector[i].end(), repptsvector[j]);     	    
            return ok;
        }
}

void MarkersCluster::assignlqmarkers(list<int> & lqmarkers, vector<list<int> > & clusters, vector<list<int> > & reppts, int * markerassignedarray, double lowqualc)
{
    cout << "number of low-quality markers ="<<lqmarkers.size()<<endl;
    int numclusters = clusters.size();
    cout << "numclusters = "<<numclusters<<endl;
    int nreppts;
    int constant=lowqualc;
    int marker;
    double lod, loddiff;
    double maxlod=std::numeric_limits<double>::min();
    int maxlodpos=-1, maxlodreppt=-1;
    double secondmaxlod=std::numeric_limits<double>::min();
    int secondmaxlodpos=-1;
    int numassigned = 0;
    int numnotassigned = 0;
    vector<list<int>> singletons; 
    list<int> singletonmarkers;

    /* make sure the number of clusters matches the number of reppt lists */
    if(numclusters != reppts.size()){
        cout << "error: number of clusters does not match number of rep. point lists"<<endl;
    }
    
    std::list<int>::iterator it = lqmarkers.begin();
    int lqcount=0;
    while(it!=lqmarkers.end()){ //iterate through all low-quality markers
        lqcount++;
	marker = *it;
        maxlodpos=-1;
        maxlod=std::numeric_limits<double>::min();
        secondmaxlodpos=-1;
        secondmaxlod=std::numeric_limits<double>::min();
        //iterate over all the clusters
        for(int cnum=0; cnum < clusters.size(); cnum++){
            //iterate over each reppt in a cluster
            for(std::list<int>::iterator repiter=(reppts[cnum]).begin(); repiter!=(reppts[cnum]).end(); repiter++){
                lod = get<2>(returnlod(marker,*(repiter)));
		if(lod>maxlod && maxlodpos != cnum){
	           maxlod = lod;
		   maxlodpos = cnum;
		   maxlodreppt = *(repiter);
		} else if (lod>secondmaxlod && cnum !=maxlodpos && secondmaxlodpos != cnum){
		   secondmaxlod = lod;
		   secondmaxlodpos = cnum;
		}
            }
        }
        
        loddiff = maxlod - secondmaxlod;
        if(loddiff > constant){
            //cout << "adding to cluster " << maxlodpos << endl;
            (clusters[maxlodpos]).push_back(marker);
            numassigned++;
	    //markerassignedarray[marker] = maxlodreppt;
        } else {
            numnotassigned++;
  	    list<int> singleton;
            singleton.push_back(marker);
            singletonmarkers.push_back(marker);
            singletons.push_back(singleton);
	}
        it++;
	if((lqcount%1000)==0) {
		cout << "processed "<< lqcount << " low-quality markers "<<endl;
                cout << "number of low-quality markers assigned = " << numassigned << endl;
    		cout << "number of low-quality markers not assigned = " << numnotassigned <<endl;
    	}
    }
    lqmarkers = singletonmarkers;
    cout << "number of low-quality markers assigned = " << numassigned << endl;
    cout << "number of low-quality markers not assigned = " << numnotassigned <<endl;
    clusters.insert(clusters.end(), singletons.begin(), singletons.end());
}

MarkersCluster::MarkersCluster(ifstream & input1, char* argv[], double sens):sensitivity(sens)
{
	string line;
	size_t maxsnps = 0;
        int ninds=0, indsunknown=1;
        
	while(getline(input1,line))
	{
                vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t "));
		markers1.push_back(strs[0]);
		calls1.push_back(vector<sparsecall>());	// push empty placeholder vector
		vector<sparsecall> mycalls;
		for(size_t i=1; i< strs.size(); ++i)
		{
			int val = atoi(strs[i].c_str());
			if(val == -1 || val == 1) // nonzero [toss 2 or 3 here]
			{
				mycalls.push_back(sparsecall(i, val));
			}
 		        else if(indsunknown && (val == -1 || val == 1 || val == 0)){
 		            	//ninds++;
				ninds = strs.size();
 			}
		}
		maxsnps = std::max(mycalls.size(), maxsnps);
                calls1.back().swap(mycalls);
	        if(indsunknown) indsunknown = 0;
	}
	
	cout << "Read input file "<< endl;
	cout << "markers.size() = " << markers1.size()<<endl;
	cout << "1st 10 marker names in input file: "<<endl;
  	for(int x = 0; x<10; x++){
          cout << markers1[x] <<endl;
	}

        cout << "Symmetric data" << endl;
	markers2 = markers1; // ABAB: Deep copy, optimize later
	calls2 = calls1;
	nmarkers = calls1.size();
	nindividuals = ninds-1;
	cout << "nindividuals = "<<nindividuals<<endl;
}

/* isBoundary(reppointlist, highest-lod-scoring-reppt, marker, maxlod-for-marker, threshold):
 *   checks whether marker should be a new representative point in repptlist;
 *   if so, adds it to the correct end; if not, does nothing
 **/
void MarkersCluster::checkBoundary(int clusternum, std::list<int> & reppts, int highestreppt, int marker, double maxlod, double threshold)
{
    int r2;
    double lodr1r2, lodr2;
    double lodmr1 = maxlod;
   
    //cout << " in checkBoundary: checking boundary points for cluster "<<clusternum<<" and marker "<<marker <<" and highestreppt "<<highestreppt << endl;
   if(highestreppt == marker){
      cout << "DUPLICATE! "<<endl;
      cout << "marker "<<marker<<" matches highest rep pt "<<highestreppt << " in cluster "<<clusternum<<endl;
      return;
   } 
   
   if(reppts.size()==1){
        //cout << "reppts.size() == 1; reppts.front() = "<<reppts.front()<<endl;
        reppts.push_back(marker);
	//cout << "added marker "<<marker<<"to back of cluster "<<clusternum<<"whose repptlist is:"<<reppts.front()<<endl;
   }else if((highestreppt == reppts.front()) || (highestreppt == reppts.back()) ){
        //determine whether the added marker should be a new representative point
        if(highestreppt == reppts.front()){
            //cout << "\tpotential new boundary point at front"<<endl;
	    //cout << "\tcluster "<<clusternum<<".front() = "<<reppts.front()<<endl;
            if(reppts.size()>1){
                std::list<int>::iterator itr = reppts.begin();
                itr++;
                r2 = *itr;
            } else {
                r2 = highestreppt;
            }
            lodr2 = get<2>(returnlod(r2,marker));
            lodr1r2 = get<2>(returnlod(r2,highestreppt));
            if( (lodr1r2 > lodr2) && (maxlod>lodr2) ){
                //cout << " adding new reppoint"<<marker <<" to beginning of cluster "<< clusternum  << endl;
                //cout << "lodr1r2 = "<<lodr1r2<<", lodr2="<<lodr2<<", maxlod = "<<maxlod<<endl;
                reppts.insert(reppts.begin(),marker);
                //cout << "added reppoint "<< m <<" to cluster " << maxlodpos << endl;
                //eliminate closest rep pt if next-to-closet one is still above threshold
                if((lodr2 > threshold) && (r2 != highestreppt)){
                    reppts.remove(highestreppt);
		    //cout <<"removed rep pt"<<highestreppt<<endl;
                }
            }
        } else {
            //cout << "\tpotential new boundary point at back"<<endl;
	    //cout << "\tcluster "<<clusternum<<".back() = "<<reppts.back()<<endl;
            std::list<int>::iterator bdryit2 = reppts.end();
            bdryit2--;
	    //cout << "*bdryit2 = "<<*bdryit2<<endl;
            bdryit2--;
	    //cout << "*bdryit2 = "<<*bdryit2<<endl;
            r2 = *bdryit2;
	    //cout << "r2 = "<<r2<<endl;;
            lodr1r2 = get<2>(returnlod(r2,highestreppt));
            lodr2 = get<2>(returnlod(r2,marker));
            if( (lodr1r2 > lodr2) && (maxlod>lodr2)){
                //cout << "adding new reppoint "<<marker <<"to end of cluster "<<clusternum<<endl;
                //cout << "lodr1r2 = "<<lodr1r2<<", lodr2="<<lodr2<<", maxlod = "<<maxlod<<endl;
                reppts.insert(reppts.end(),marker);
                //eliminate closest rep pt if next-to-closest one is still above threshold				   
                if(lodr2>threshold){				   
                    reppts.remove(highestreppt);
		    //cout <<"removed rep pt"<<highestreppt<<endl;
                }
            }
        } 
    }
}


/* bubbleCluster: clusters a set of markers by
 *  the BubbleCluster algorithm
 **/
void MarkersCluster::bubbleCluster(vector<list<int>> & clustervec, vector<list<int>> & repptvec, int *markerarray,int nmarkers,
                   int *markerassigned, double threshold, double selflodthreshold)
{
    
    int numassignedmarkers=0;
    int nclusters  = 0;
    int nclusterscreated = 0;
    
    //Set the first cluster to contain the first random marker (in range 0...#markers-1)
    // and the first representative point list to also contain just the first random marker
    std::list<int> cluster0;
    cluster0.push_back(*markerarray);
    clustervec[0] = cluster0;
    repptvec[0] = cluster0;
    nclusters++;
    nclusterscreated++;
    
    double lod, maxlod, minlod=std::numeric_limits<double>::max();
    double secondmaxlod = std::numeric_limits<double>::min(),prevmaxlod,lodr1r2, lodr2;
    int m,j,mcounter=1, nummaxlodscounter=0, nummatchingclusterassignments=0, selflodok=1;
    int maxlodpos=-1, secondmaxlodpos=-1, maxlodreppoint, r2;
    pair<int,int> initpair(-1,-1);
    pair<int,int> newpair(-1,-1);
    vector<int> matchingclusternumbers(20,-1);
    int nummatching=0; 
    double nummatchingavg=0;
	
    cout << "starting clustering process..." << endl;
    while((mcounter < nmarkers) && selflodok){
        //pick a marker
        m = *(markerarray+mcounter);
        mcounter++;
        // if the marker is of low quality, then all subsequent ones will be too;
        // break and place the rest of the markers according to new rules
        if(get<2>(returnlod(m,m)) <= selflodthreshold){
            cout << "self-LOD less than "<<selflodthreshold<<" at "<<mcounter<< " markers; breaking loop..."<<endl;
            selflodok=0;
	    break;
        }
        
        maxlod=std::numeric_limits<double>::min();
        nummaxlodscounter=0;
        maxlodpos=-1;
        secondmaxlodpos=-1;
        secondmaxlod=std::numeric_limits<double>::min();
        nummatching=0;
        for(int l=0; l<matchingclusternumbers.size(); l++){    
	  matchingclusternumbers[l] = -1;	
        }
	
        for(int i=0; i<repptvec.size(); i++){
            //for each cluster, find the representative point that has a highest lod to the marker
            // also find all representative points for which the LOD score is higher than the threshold
	    for(std::list<int>::iterator repptiter= repptvec[i].begin();repptiter!=repptvec[i].end(); repptiter++){
                lod = get<2>(returnlod(m,*(repptiter)));
	        
                if(lod>maxlod){
		  if(lod>threshold){
                    if(nummatching==0){
                        matchingclusternumbers[nummatching] = i;
                        nummatching++;
	            } else if(matchingclusternumbers[(nummatching-1)] != i){
                        matchingclusternumbers[nummatching] = i;
                        nummatching++;
                    } 
                  }
		  maxlod = lod;
                  maxlodpos=i;
                  maxlodreppoint = *(repptiter);
                } else if(lod>threshold && maxlodpos!=i && maxlodpos!=-1){
		// in this case, lod>threshold but lod <=maxlod
                   if(nummatching==0){
                       matchingclusternumbers[nummatching] = i;
		       nummatching++;
		    } else if(matchingclusternumbers[(nummatching-1)] != i){
		      matchingclusternumbers[nummatching] = i;
		      nummatching++;
		    }
		    if(lod>secondmaxlod && lod <= maxlod && secondmaxlodpos != i && maxlodpos != i){
                        secondmaxlod=lod;
                        secondmaxlodpos=i;
                    }
		} 
	   }
        }
        // the number of the cluster with max lod to the marker must be >= 0 
        assert(maxlodpos>=0);
        
        //if the maximum lod score is less than a threshold,
        //create a new cluster with this marker as the sole point in the cluster
        //and the sole point in the representative point list for this cluster
	nummatchingavg += (double) nummatching;
	int l=0;
        if(maxlod < threshold){
            std::list<int> cluster;
            cluster.push_back(m);
            clustervec[nclusterscreated] = cluster;
            repptvec[nclusterscreated] = cluster;
            nclusters++;
            nclusterscreated++;
           // *(markerassigned + m) = m;
            //cout << " created new cluster with marker "<< m <<"; nclusters = "<<nclusters<<endl;
            //cout << " nclusterscreated = "<<nclusterscreated<<endl;
        }
        
	// otherwise, add the marker to the cluster with max lod score
        // and check if it is a new boundary point
        else{
	  //the maxlod between the marker and at least one rep pt is higher than the threshold
	  // add the marker to the cluster with the rp pt that scored highest with the marker
	  // and check the boundary (check whether m should be a new boundary pt)
          (clustervec[maxlodpos]).push_back(m);
          checkBoundary(maxlodpos,repptvec[maxlodpos],maxlodreppoint, m, maxlod, threshold);
            
        //if the marker has a lod higher than threshold to more than one cluster, merge the two clusters
          if(nummatching>1){
            sort(matchingclusternumbers.begin(), matchingclusternumbers.end(),greaterthan);
            unique(matchingclusternumbers.begin(), matchingclusternumbers.end());
	    l=0;
	    while(matchingclusternumbers[l]>0){
              if(matchingclusternumbers[l] != maxlodpos){
                if( 0 == mergerepptlists(repptvec,maxlodpos,matchingclusternumbers[l],threshold, maxlod) ){
                  clustervec[maxlodpos].splice(clustervec[maxlodpos].end(),clustervec[matchingclusternumbers[l]]);
                  //cout << "merged clusters "<<maxlodpos<<" and "<< matchingclusters[l].first <<endl;
                  nclusters--;
                } 
 	      }
	      l++;
	    }
	}
       }
       if(mcounter%1000 == 0){
            cout << "processed "<<mcounter<<" markers, number linkage groups = "<<nclusters<<endl;
	    cout << "number of clusters created = "<<nclusterscreated<<endl;
            cout << "self-LOD of current marker: " << get<2>(returnlod(m,m)) <<endl;
       }
    }//end while
    nummatchingavg = nummatchingavg/mcounter;
    cout << "number of matches on average: " << nummatchingavg << endl;		
}




/* mergesmallclusters:
 *    attempts to merge all small clusters with larger clusters
 **/
void MarkersCluster::mergesmallclusters(vector<list<int>> & clusters, vector<list<int>> & repptlists, double sigma, double lowqualc, double threshold)
{
    int nmarkersincluster, m1,m2,marker, maxclust = -1;
    int nummergesreppts=0, reverse=0;
    int maxreppti=-1, maxrepptj=-1, iter=0;
    int frontiter=1;
    double lodprev, loddif;
    vector<list<int>> smallclusters;
    vector<list<int>> smallclusterreppts;
    vector<int> csizes;
    int numsmallcs = 0, smallcsize;
    
    /* store all cluster sizes */
    for(int ccounter = 0; ccounter<clusters.size(); ccounter++ ){
       if(clusters[ccounter].size()>0){
         csizes.push_back(clusters[ccounter].size());
       }
    }
    /* sort the cluster sizes from greatest to smallest */
    sort(csizes.begin(), csizes.end(),greaterThan);
    
    /* determine the cut-off for a small cluster;
     * i.e. find smallcsize such that any cluster
     * with size less than smallcsize is added to
     * the small cluster list
     * */
    int firstscsize=csizes.front();
    //cout << "firstscsize = " << firstscsize << endl;
    int prevsize = firstscsize;
    int currentsize, smallcsizeindex;
    for(int i=1; (i<csizes.size() && csizes[i]>0); i++){
        currentsize=csizes[i];
        if(currentsize < (prevsize/2)){
            smallcsize = currentsize;
            smallcsizeindex = i;
            break;
	}
    }
    cout << "smallcsize found by heuristics = " <<smallcsize<<endl;
    smallcsize = sigma;
    cout << "smallcsize set to: " <<sigma<<endl;
    cout << "number of clusters with size greater than "<<smallcsize<<" = "<< (smallcsizeindex) <<endl;
    
    /* add all clusters with size less than
     * smallcsize to the small cluster list */
    for(int ccounter=0, rcounter = 0; ccounter<clusters.size(), rcounter<repptlists.size(); ccounter++, rcounter++ ){
        if(clusters[ccounter].size() <= smallcsize && clusters[ccounter].size()>0){
            list<int> sc = clusters[ccounter];
            smallclusters.push_back(sc);
  	    cout << "added small cluster of size "<<sc.size()<<" to small clusters"<<endl;
            list<int> scrp = repptlists[rcounter];
            smallclusterreppts.push_back(scrp);
            numsmallcs++;
        }
    }
    
    /* remove all clusters with size less than smallcsize
     *  from the main cluster list  */
    for(int ccounter=0, rcounter=0; ccounter<clusters.size(), rcounter<repptlists.size(); ccounter++, rcounter++){
        if(clusters[ccounter].size() <= smallcsize){
            clusters[ccounter].clear();
            repptlists[ccounter].clear();
        }
    }
    
    cout << "number of small clusters = "<<smallclusters.size()<<endl;
    cout << "number of large clusters = "<<clusters.size()<<endl;
    
    /* loop through all reppts in all large clusters and compare them to
     * a random marker in each small cluster. If the highest LOD score between
     * the random marker and a rep pt. in cluster A is greater than a
     * constant+the highest LOD score to a rep pt. in any other cluster B,
     * then merge the small cluster with cluster A*/
    
    double lodsc, lodscfront,lodscback,maxlod, secondmaxlod;
    int i=0, j=0;
    int maxlodpos=-1;
    int secondmaxlodpos=-1;
    
    for(int sccount=0; sccount<smallclusters.size(); sccount++){
        maxlod=numeric_limits<double>::min();
        secondmaxlod=numeric_limits<double>::min();
        maxlodpos=-1;
	secondmaxlodpos=-1; 
	i=0;
        for(int ccounter=0, rcounter=0; ccounter < clusters.size(), rcounter<repptlists.size(); ccounter++, rcounter++){
            if(clusters[ccounter].size()>0 && repptlists[ccounter].size()>0 && repptlists[ccounter].front() != -1){
                for(list<int>::iterator rpit = repptlists[ccounter].begin(); rpit != repptlists[ccounter].end(); rpit++){
                    //select 5 random points from the small cluster if its size is >= 5
                    if(smallclusters[sccount].size() >= 5){
                      double maxlodsc=0;
                      for(int si=0;si<5; si++){
                       int randomMarkerIndex = rand()%smallclusters[sccount].size();
                       auto rmit=smallclusters[sccount].begin();
		       int rmi=0;
                       while(rmi<randomMarkerIndex){
			 rmit++;
		       	 rmi++;
		       }
                       //cout << "random marker selected from small cluster = "<<*rmit<<endl;
		       double lod = get<2>(returnlod(*rmit,*rpit));
		       if(lod>maxlodsc){ 
			 maxlodsc = lod;
		       }
 		      }
                      lodsc = maxlodsc;
		    } else { //otherwise, use all points in the small cluster
                      double maxlodsc=0;
                      for(auto smallitr=smallclusters[sccount].begin(); smallitr!=smallclusters[sccount].end(); smallitr++){
		         double lod = get<2>(returnlod(*smallitr,*rpit));
		         if(lod>maxlodsc) maxlodsc = lod;
                      }
                      lodsc = maxlodsc;
		    }
		    if(maxlodpos!=i && lodsc>maxlod ){
                            maxlod=lodsc;
                            maxlodpos=i;
                    }else if(lodsc>secondmaxlod && secondmaxlodpos != i && maxlodpos != i){
                        secondmaxlod = lodsc;
                        secondmaxlodpos=i;
                    }
                }
            }
            i++;
        }
        double c = lowqualc;
        if((maxlod - secondmaxlod) > c && (maxlodpos != j) && (maxlodpos != -1)){
            clusters[maxlodpos].insert(clusters[maxlodpos].end(),(smallclusters[sccount]).begin(), (smallclusters[sccount]).end());
	    repptlists[maxlodpos].insert(repptlists[maxlodpos].end(), (smallclusterreppts[sccount]).begin(), (smallclusterreppts[sccount]).end());
            smallclusters.erase(smallclusters.begin()+sccount);
            smallclusterreppts.erase(smallclusterreppts.begin()+sccount);
	    sccount--;
	} else {
            //cout << "could not place small cluster " << j << endl;
            //cout << "\tmaxlod = " <<maxlod<<", secondmaxlod = "<<secondmaxlod<< endl;
        }
        j++; 
    }
    clusters.insert(clusters.end(),smallclusters.begin(), smallclusters.end());
    repptlists.insert(repptlists.end(), smallclusterreppts.begin(), smallclusterreppts.end());
}


/* printclusters(clusterlist, outputfilename, vectorofmarkernames):
 *    prints out one file per cluster found, with one line per marker in that cluster
 **/
void MarkersCluster::printclusters(vector<list<int>> & clusters, char* ofname)
{    
    int nclusters = clusters.size();
    stringstream ss;
    string cnumstr;
    int numclusterswith1marker=0;
    int numclusterswith2markers=0;
    int numclusterswithmorethan10markers=0;
    int numclusterswithmorethan100markers=0;
    int numclusterswithmorethan1000markers=0;
    int numprintedclusters=0;
    vector<int> singletons;
    for(int i=0; i<nclusters; i++)
    {
        if(clusters[i].size()>1) {
            ss << ofname << "v4cluster" << numprintedclusters;
            cnumstr=ss.str();
            if((clusters[i]).size() == 1) numclusterswith1marker++;
            else if((clusters[i]).size() == 2) numclusterswith2markers++;
            else if((clusters[i]).size() > 1000) numclusterswithmorethan1000markers++;
            else if((clusters[i]).size() > 100) numclusterswithmorethan100markers++;
            else if((clusters[i]).size() > 10) numclusterswithmorethan10markers++;
            ofstream output(cnumstr);
            list<int>::iterator itout = clusters[i].begin();
            while(itout!=clusters[i].end())
            {
                output << markers1[*itout] << endl;
                itout++;
            }
            output.close();
            numprintedclusters++;
        } else if (clusters[i].size() == 1){
	    singletons.push_back(clusters[i].front());
	}
        ss.str("");
    }
   
    ss.str("");
    ss << "singletons_"<<ofname << "v4";
    cnumstr=ss.str();
    ofstream outputsingletons(cnumstr);
    for(int j = 0; j<singletons.size(); j++){
	outputsingletons << markers1[singletons[j]] << endl;
    } 
    outputsingletons.close();

    printf("number of clusters with 1 marker = %d\n",numclusterswith1marker);
    printf("number of clusters with 2 markers = %d\n",numclusterswith2markers);
    printf("number of clusters with 10-99 markers = %d\n", numclusterswithmorethan10markers);
    printf("number of clusters with 100-999 markers = %d\n", numclusterswithmorethan100markers);
    printf("number of clusters with more than 1000 markers = %d\n",numclusterswithmorethan1000markers);
}

/* printreppoints(repptlist, outputfilename, vectorofmarkernames):
 *    prints out one file per cluster, with one line per representative point for that cluster
 **/
void MarkersCluster::printreppoints(vector<list<int>> & reppointvec, char* ofname)
{
    stringstream ssrp;
    string rpnumstr;
    int numprintedreppts=0;
    
    for(int i=0; i<reppointvec.size(); i++)
    {
        std::list<int> rpl = reppointvec[i];
        if(rpl.size()>0)
        {
            ssrp << ofname << "v4reppoints" << i;
            rpnumstr = ssrp.str();
            ofstream outputrp(rpnumstr);
            for(std::list<int>::iterator itout=rpl.begin(); itout!=rpl.end(); ++itout){
                outputrp << markers1[*itout] << endl;
            }
            outputrp.close();
            ssrp.str("");
            numprintedreppts++;
        }
        else
        {
            ssrp.str("");
        }
    }
    cout << "number of printed reppoint lists = "<<numprintedreppts<<endl;
}

void MarkersCluster::printMarkerRptAssignments(int *markerassignedarray, int nummarkers){
   ofstream outf("markerassignments.out");
   outf << "MarkerNumber\tRepPtAssignment"<<endl;
   for(int i=0; i<nummarkers; i++){
      outf << markers1[i] <<"\t"  << markers1[markerassignedarray[i]] <<endl;
   }
   outf.close();
   cout << "printed marker assignments to rep points."<<endl;
}

/* printRptLODMtx: print out a matrix of LOD scores between all
 * representative (sketch) points in cluster number cnum
 */
void MarkersCluster::printRptLODMtx(list<int> & rpptlist, int cnum){
   stringstream ss;
   ss << "rpptlist"<<cnum<<".mtx"; 
   ofstream outf(ss.str());
   for(list<int>::iterator it = rpptlist.begin(); it!=rpptlist.end(); it++){
    for(list<int>::iterator it2 = rpptlist.begin(); it2!=rpptlist.end(); it2++){
       outf <<get<2>(returnlod(*(it),*(it2)))<<" ";	
    }
    outf<<endl;
   }
   outf.close();
}

/* getNumNonmissing: return the number of non-zero entries for a
 * particular marker referenced by marker number markernum. I.e.,
 * this function returns the number of non-missing values for a 
 * particular marker. */
int MarkersCluster::getNumNonmissing(int markernum){
   vector<sparsecall> markervec = calls1[markernum];
   return markervec.size();
}
