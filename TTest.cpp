//
//  TTest.cpp
//  
//
//  Created by Kim Ngo Aug 11, 2014.

// 	Calculates TTest for regions of a vector of floats
//	Confidence level 95%
//	http://en.wikipedia.org/wiki/Welch%27s_t_test

#include "TTest.h"
#include <vector>
#include <math.h>
#include <cmath>		// abs
#include <utility> 	// pair
#include <iostream> 

// Calculate T-Test
TTest :: TTest(const vector<float>& mm_seq, const pair<int, int>& region1, const pair<int, int>& region2) {

//	cout << "Region A: " << region1.first << "," << region1.second << endl;
//	cout << "Region B: " << region2.first << "," << region2.second << endl;

	region1_len = region1.second - region1.first +1;
	region2_len = region2.second - region2.first +1;
	
//	cout << "N\t" << region1_len << " " << region2_len << endl;
	region1_avg = calc_region_mm_avg(mm_seq, region1.first, region2.second);
	region2_avg = calc_region_mm_avg(mm_seq, region2.first, region2.second);
	
//	cout << "avg\t" << region1_avg << " " << region2_avg << endl;
	region1_var = calc_region_var(mm_seq, region1.first, region1.second, region1_avg);
	region2_var = calc_region_var(mm_seq, region2.first, region2.second, region2_avg);
	
//	cout << "var\t" << region1_var << " " << region2_var << endl;
	df = calc_df(region1_var, region1_len, region2_var, region2_len);
	t = calc_t(region1_avg, region1_var, region1_len, region2_avg, region2_var, region2_len);
	sig = calc_significance(t, df);
}

float TTest :: get_avg1() {
	return region1_avg;
}
float TTest :: get_avg2() {
	return region2_avg;
}

// Returns greater average to determine slow down or increase
int TTest :: greater() {
	if (get_avg1() > get_avg2()) return 1;
	else return 2;
}


// Calculates average %Min-Max value for region
float TTest :: calc_region_mm_avg(const vector<float>& mm_seq, const int& begin, const int& end) {

	float nan = 0;
	float sum = 0;
	for (int i = begin; i < end; i++) {
		if (isnan(mm_seq[i])) {
			nan++;
			continue;
		}
		sum += mm_seq[i];
	}
	float n = end - begin - nan;
	return sum/n;
}
// Calculates standard deviation for %Min-Max value for region
float TTest :: calc_region_var(const vector<float>& mm_seq, const int& begin, const int& end, const float& avg) {
	
	float nan = 0;
	float var = 0;
	for (int i = begin; i < end; i++) {
		if (isnan(mm_seq[i])) {
			nan++;
			continue;
		}
		var = pow((mm_seq[i] - avg), 2);
	}
	float n = end - begin - nan;	// n number of values
	return var/(n-1);
}
// Calculates degrees of freedom
float TTest :: calc_df(const float& var1, const float& n1, const float& var2, const float& n2) {
	float numerator = pow( var1/n1 + var2/n2, 2);
	float denominator = pow(var1,2)/(pow(n1,2)*(n1-1)) + pow(var2,2)/(pow(n2,2)*(n2-1));
//	cout << "DF = " << numerator/denominator << endl;
	return numerator/denominator;
}
// Calculates statistic t
float TTest :: calc_t(const float& avg1, const float& var1, const float& n1, const float& avg2, const float& var2, const float& n2) {
	float numerator = abs(avg1 - avg2);
	float denominator = pow( (var1/n1 + var2/n2), 0.5);
//	cout << avg1 << "," << var1 << "," << n1 << endl;
//	cout << avg2 << "," << var2 << "," << n2 << endl;
	cout << "T = " << numerator << "/" << denominator << " = " << numerator/denominator << endl;
	return numerator/denominator;
}


bool TTest :: significant() {
	return sig;
}

// Returns true of significantly different
bool TTest :: calc_significance(float& t, float& df) {

	int degrees = ceil(df);
	if (degrees > 30) {
		if (degrees > 30 && degrees <= 40) {
			degrees = 40;
		} else if (degrees > 40 && degrees <= 60) {
			degrees = 60;
		} else if (degrees > 60 && degrees <= 80) {
			degrees = 80;
		} else if (degrees > 80 && degrees <= 100) {
			degrees = 100;
		} else degrees = 1000;
	}
	if (isnan(df)) return false;
//	cout << df << "->" << degrees << endl;

	switch(degrees) {
		case 1:
			if (t > 12.71) return true;
		case 2:
			if (t > 4.303) return true;
		case 3:
			if (t > 3.182) return true;
		case 4:
			if (t > 2.776) return true;
		case 5:
			if (t > 2.571) return true;
		case 6:
			if (t > 2.447) return true;
		case 7:
			if (t > 2.365) return true;
		case 8:
			if (t > 2.306) return true;
		case 9:
			if (t > 2.262) return true;
		case 10:
			if (t > 2.228) return true;
		case 11:
			if (t > 2.201) return true;
		case 12:
			if (t > 2.179) return true;
		case 13:
			if (t > 2.160) return true;
		case 14:
			if (t > 2.145) return true;
		case 15:
			if (t > 2.131) return true;
		case 16:
			if (t > 2.120) return true;
		case 17:
			if (t > 2.110) return true;
		case 18:
			if (t > 2.101) return true;
		case 19:
			if (t > 2.093) return true;
		case 20:
			if (t > 2.086) return true;
		case 21:
			if (t > 2.080) return true;
		case 22:
			if (t > 2.074) return true;
		case 23:
			if (t > 2.069) return true;
		case 24:
			if (t > 2.064) return true;
		case 25:
			if (t > 2.060) return true;
		case 26:
			if (t > 2.056) return true;
		case 27:
			if (t > 2.052) return true;
		case 28:
			if (t > 2.048) return true;
		case 29:
			if (t > 2.045) return true;
		case 30:
			if (t > 2.042) return true;
		case 40:
			if (t > 2.021) return true;
		case 60:
			if (t > 2.000) return true;
		case 80:
			if (t > 1.990) return true;
		case 100:
			if (t > 1.984) return true;
		case 1000:
			if (t > 1.962) return true;
		default:
			return false;
	}
}
