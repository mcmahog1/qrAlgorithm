
#include <Rcpp.h>
#include<vector>
#include <math.h> 
using std::vector ;

using namespace Rcpp ;



//[[Rcpp::export]]
double getInnerProduct(std::vector<double> x,std::vector<double> y) {

int x_length = x.size() ;
int y_length = y.size() ;
double product = 0 ;

	for(int i = 0 ; i < x_length ; i++) {
	
	product = product + x[i]*y[i] ;

	}


return(product) ;

}


//[[Rcpp::export]]
double getLength(std::vector<double> x) {

int x_length = x.size() ;
double length = 0 ; 
	for(int i = 0 ; i < x_length ; i++) {

	length = length + x[i]*x[i] ;

	}

length = sqrt(length) ;

return(length) ;

}



//[[Rcpp::export]]
std::vector<double> getProjection(std::vector<double> x, std::vector<double> y) {
int y_length = y.size() ;

double numerator = getInnerProduct(x,y) ;
double denominator = getInnerProduct(y,y) ;

double ratio = numerator/denominator ;

	for(int i = 0 ; i < y_length ; i++) {

	y[i] = y[i]*ratio ;

	}

return(y) ;
}



//[[Rcpp::export]]
NumericMatrix getU(NumericMatrix x) {

int rows = x.nrow() ;
int cols = x.ncol() ;

NumericMatrix U(rows,cols) ;

for(int i = 0 ; i < cols ; i++ ) {

	NumericVector myU = x( _, i) ;
	std::vector<double> myU_vec = Rcpp::as<std::vector<double> >(myU) ;

		for(int j = 0 ; j < i ; j++ ) {
	
			NumericVector a = x( _, i);
			std::vector<double> a_vec = Rcpp::as<std::vector<double> >(a) ;

			NumericVector e = U( _, j);
			std::vector<double> e_vec = Rcpp::as<std::vector<double> >(e) ;

			std::vector<double> projected = getProjection(a_vec,e_vec) ;
			
				for(int m = 0 ; m < rows ; m++) {				

				myU[m] = myU[m] - projected[m] ;

				}

		}

	for(int k = 0 ; k < rows ; k++) {

	U(k,i) = myU[k] ;

	}	


}


return(U) ;

}


//[[Rcpp::export]]
NumericMatrix getQ(NumericMatrix x) {

int rows = x.nrow() ;
int cols = x.ncol() ;

	for(int i = 0 ; i < cols ; i++) {

	NumericVector col = x(_,i) ;
	std::vector<double> col_vec = Rcpp::as<std::vector<double> >(col) ;
	double length = getLength(col_vec) ;
	
		for(int j = 0 ; j < rows ; j++ ) {

		x(j,i) = x(j,i)/length ;

		}

	}

return(x) ;
}


//[[Rcpp::export]]
NumericMatrix getR(NumericMatrix Q, NumericMatrix A) {

int rowsQ = Q.nrow() ;
int colsQ = Q.ncol() ;

int rowsA = A.nrow() ;
int colsA = A.ncol() ;

NumericMatrix R(colsQ,colsA) ;

std::fill(R.begin(),R.end(),0) ;

	for(int i = 0 ; i < colsQ ; i++) {

	NumericVector colQ = Q(_,i) ;
	std::vector<double> colQ_vec = Rcpp::as<std::vector<double> >(colQ) ;
	
	
		

		for(int j = 0 ; j < colsA ; j++) {

		NumericVector colA = A(_,j) ;
		std::vector<double> colA_vec = Rcpp::as<std::vector<double> >(colA) ;
	
		double product = getInnerProduct(colA_vec,colQ_vec) ;

		R(i,j) = product ;

		}


	}

return(R) ;

}
