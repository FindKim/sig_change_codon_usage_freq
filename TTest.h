//
//  TTest.h
//  
//
//  Created by Kim Ngo Aug 11, 2014.

// 	Calculates TTest for regions of a vector of floats
//	Confidence level 95%
//	http://en.wikipedia.org/wiki/Welch%27s_t_test

#ifndef TTEST_H
#define TTEST_H

#include <string>
#include <vector>
#include <utility>

using namespace std;

class TTest {

	public:
		TTest(const vector<float>&, const pair<int, int>&, const pair<int, int>&);
			// vector of mm values
			// begin & end for region 1 & region 2
		bool significant();
			// Returns significance
		int greater();
			// Returns which average is greater
			// 1 for first, 2 for seond
			
	private:
		float region1_len;
		float region2_len;
		float region1_avg;
		float region2_avg;
		float region1_var;
		float region2_var;
		float df;
		float t;
		bool sig;
		
		float get_avg1();
		float get_avg2();
		bool calc_significance(float&, float&);
			// t & df
			// Returns true if significant difference with confidence of 95%
		float calc_t(const float&, const float&, const float&, const float&, const float&, const float&);
			// avg1, var1, n1, avg2, var2, n2
			// Returns t score
		float calc_df(const float&, const float&, const float&, const float&);
			// var1, n1, var2, n2
			// Returns degrees of freedom
		float calc_region_var(const vector<float>&, const int&, const int&, const float&);
			// mm_seq, begin, end, avg
			// Returns mm variance for region
		float calc_region_mm_avg(const vector<float>&, const int&, const int&);
			// mm_seq, begin, end
			// Returns mm avgerage for region
		
};

#endif

