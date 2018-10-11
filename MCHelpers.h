/* Copyright 2014 Aydin Buluc, Veronika Strnadova-Neeley
 * 
 *  MCHelpers.h is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MCHelpers.h is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MCHelpers.h.  If not, see <http://www.gnu.org/licenses/>.  
 */
#ifndef _MC_HELPERS_
#define _MC_HELPERS_

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

using namespace std;

struct sparsecall
{
    sparsecall(): laneid(0), phase(0){};
	sparsecall(int mylane, int myphase): laneid(mylane), phase(myphase) {};
	
	int laneid;
	int phase; // -1 or 1
	int numagree;
        int numdisagree;
        	
	bool operator== (const sparsecall & rhs)
	{
		return (laneid == rhs.laneid);
	}
	bool operator< (const sparsecall & rhs)
	{
		return (laneid < rhs.laneid);
	}
	bool operator> (const sparsecall & rhs)
	{
		return (laneid > rhs.laneid);
	}
};

ostream & operator<< ( ostream & out, const sparsecall & spcall )
{
    out << "(" << spcall.laneid << "," << spcall.phase << ")";
    return out;
};
        
// only called on items with same laneid's (for merging)
struct multcalls
{
    int operator()(const sparsecall & x, const sparsecall & y) const
    {
        assert( (x.laneid == y.laneid) );
        return x.phase * y.phase;  // returns 1 if same phase, -1 if opposite phase
    }
};

// returns of pair of <flops, accum> where
// - flops is the size of the intersection set
// - accum is the accumulated value as defined by __binop
// (binop here counts each same sign as 1 and each different sign as -1)
template<typename _InputIterator1, typename _InputIterator2, typename _BinOp>
pair<int,int>  intersect_accumulate(_InputIterator1 __first1, _InputIterator1 __last1,
                                    _InputIterator2 __first2, _InputIterator2 __last2, _BinOp __binop)
{
    typedef typename iterator_traits<_InputIterator1>::value_type _ValueType1;
    typedef typename iterator_traits<_InputIterator2>::value_type _ValueType2;
            
    int accum = 0;
    int flops = 0;
    while (__first1 != __last1 && __first2 != __last2)
    {
        if (*__first1 < *__first2)
            ++__first1;
        else if (*__first2 < *__first1)
            ++__first2;
        else
        {
            ++flops;
            accum += __binop(*__first1, *__first2);
            ++__first1;
            ++__first2;
        }
    }
    return make_pair(flops, accum);
};
        
/* flip_phase: flips -1's to 1's and 1's to -1's in a 
 *   given vector of sparsecalls
 */
template<typename _InputIterator1>
void flip_phase(_InputIterator1 __first1, _InputIterator1 __last1){
    //typedef typename iterator_traits<_InputIterator1>::value_type _ValueType1;
    while (__first1 != __last1)// && __first2 != __last2)
    {
        (*__first1).phase *= -1;
	++__first1;
    }
};

bool greaterThan (int i,int j) { return (i>j); }

#endif
