/*
 * CalculateAverages
 *
 * Author:  J. Doe
 *
 * Creation Date:    1/20/95
 *	modified:    1/25/95	...(reason)...
 *
 * This function takes an array of real numbers and calculates the
 *	mean, median, and mode of those values.
 *
 * Signature:
 *	void calculateAvgs (double valArray[], int n, double & mean,
 *	  		    double & median, double & mode)
 *	pre-condition: valArray contains n values
 *	post-condition:
 *	  	mean = (sum of all n values in valArray)/n;
 *	  	mode = the number that appears most frequently ...
 *	  	if n is odd, then median is ...
 *	  	if n is even, then median is ...
 *
 * Input:
 *	The function has two input parameters: valArray, and n.
 *	The first is an array of double, which contains the
 *	values for which the function should calculate the mean,
 *	mode, and median.  The second parameter indicates how
 *	many values there are in the array.
 *
 * Output:
 *	The program calculates the mean, median, and mode of
 *	the real numbers provided as input, and passes them
 *	back to the calling function as output parameters.
 *
 * Implementation:
 *	The function first sorts the real numbers provided as
 *	input, to make it easier to calculate the median.  It
 *	then calculates the mean, median, and mode according to
 *	the following algorithms:
 *	  	Mean: (sum of values)/(number of values)
 *	  	Median: ...
 *	  	Mode: ...
 *
 */
