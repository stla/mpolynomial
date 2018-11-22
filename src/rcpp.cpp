
#include <Rcpp.h>
#include "polynomial.hpp"
using namespace std;

Rcpp::IntegerVector burkardt_unrank_grlex(int m, int rank)
{
  int i;
  int j;
  int ks;
  int nksub;
  int nm;
  int ns;
  int r;
  int rank1;
  int rank2;
  Rcpp::IntegerVector x(m);
  int *xs;
  //
  //  Ensure that 1 <= M.
  //
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
  //
  //  Ensure that 1 <= RANK.
  //
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    exit ( 1 );
  }
  //
  //  Special case M == 1.
  //
  if ( m == 1 )
  {
    x[0] = rank - 1;
    return x;
  }
  //
  //  Determine the appropriate value of NM.
  //  Do this by adding up the number of compositions of sum 0, 1, 2, 
  //  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
  //  gives you the rank of the composition within the set of compositions
  //  of sum NM.  And that's the number you need in order to do the
  //  unranking.
  //
  rank1 = 1;
  nm = -1;
  for ( ; ; )
  {
    nm = nm + 1;
    r = i4_choose ( nm + m - 1, nm );
    if ( rank < rank1 + r )
    {
      break;
    }
    rank1 = rank1 + r;
  }
  
  rank2 = rank - rank1;
  //
  //  Convert to KSUBSET format.
  //  Apology: an unranking algorithm was available for KSUBSETS,
  //  but not immediately for compositions.  One day we will come back
  //  and simplify all this.
  //
  ks = m - 1;
  ns = nm + m - 1;
  xs = new int[ks];
  
  nksub = i4_choose ( ns, ks );
  
  j = 1;
  
  for ( i = 1; i <= ks; i++ )
  {
    r = i4_choose ( ns - j, ks - i );
    
    while ( r <= rank2 && 0 < r )
    {
      rank2 = rank2 - r;
      j = j + 1;
      r = i4_choose ( ns - j, ks - i );
    }
    xs[i-1] = j;
    j = j + 1;
  }
  //
  //  Convert from KSUBSET format to COMP format.
  //
  x[0] = xs[0] - 1;
  for ( i = 2; i < m; i++ )
  {
    x[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  x[m-1] = ns - xs[ks-1];
  
  delete [] xs;
  
  return x;
}
//****************************************************************************80

int burkardt_rank_grlex ( int m, Rcpp::IntegerVector x )
{
  int i;
  int j;
  int ks;
  int n;
  int nm;
  int ns;
  int rank;
  int tim1;
  int *xs;
  //
  //  Ensure that 1 <= M.
  //
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_RANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
  //
  //  Ensure that 0 <= X(I).
  //
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_RANK_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
  //
  //  NM = sum ( X )
  //
  nm = Rcpp::sum ( x );
  //
  //  Convert to KSUBSET format.
  //
  ns = nm + m - 1;
  ks = m - 1;
  xs = new int[ks];
  xs[0] = x[0] + 1;
  for ( i = 2; i < m; i++ )
  {
    xs[i-1] = xs[i-2] + x[i-1] + 1;
  }
  //
  //  Compute the rank.
  //
  rank = 1;
  
  for ( i = 1; i <= ks; i++ )
  {
    if ( i == 1 )
    {
      tim1 = 0;
    }
    else
    {
      tim1 = xs[i-2];
    }
    
    if ( tim1 + 1 <= xs[i-1] - 1 )
    {
      for ( j = tim1 + 1; j <= xs[i-1] - 1; j++ )
      {
        rank = rank + i4_choose ( ns - j, ks - i );
      }
    }
  }
  
  for ( n = 0; n < nm; n++ )
  {
    rank = rank + i4_choose ( n + m - 1, n );
  }
  
  delete [] xs;
  
  return rank;
}
//****************************************************************************80

Rcpp::List burkardt_compress(int o1, double c1[], int e1[], int m)
{
  int get;
  int put;
  const double r8_epsilon_sqrt = 0.1490116119384766E-07;
  //
  //  Add coefficients associated with the same exponent.
  //
  get = 0;
  put = 0;
  std::vector<double> c2;
  std::vector<int> e2;
  
  while ( get < o1 )
  {
    get = get + 1;
    
    if ( 0 == put )
    {
      put = put + 1;
      c2.push_back(c1[get-1]);
      e2.push_back(e1[get-1]);
    }
    else
    {
      if ( e2[put-1] == e1[get-1] )
      {
        c2[put-1] = c2[put-1] + c1[get-1];
      }
      else
      {
        put = put + 1;
        c2.push_back(c1[get-1]);
        e2.push_back(e1[get-1]);
      }
    }
  }
  
  int o2 = put;
  //
  //  Clear out zeros and tiny coefficients.
  //
  get = 0;
  put = 0;
  
  while ( get < o2 )
  {
    if ( r8_epsilon_sqrt < fabs ( c2[get] ) )
    {
      c2[put] = c2[get];
      e2[put] = e2[get];
      put = put + 1;
    }
    get = get + 1;
  }
  
  o2 = put;
  
  Rcpp::NumericMatrix Powers(o2, m);
  for(int i=0; i<o2; i++){
    Rcpp::NumericMatrix::Row powers = Powers(i, Rcpp::_);
    powers = burkardt_unrank_grlex(m, e2[i]);
  }
  
  std::vector<double>(c2.begin(), c2.begin() + o2).swap(c2);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("coefficients") = c2, 
    Rcpp::Named("powers") = Powers);
  return out;
}
//****************************************************************************80


// [[Rcpp::export]]
Rcpp::List burkardt_polynomial_dif (
    Rcpp::NumericVector c1, Rcpp::IntegerMatrix Powers, Rcpp::IntegerVector dif)
{
  int i;
  int j;
  int m = Powers.ncol();
  int o1 = c1.size();
  int o2 = o1;
  double c2[o2];
  for ( j = 0; j < o1; j++ )
  {
    c2[j] = c1[j];
  }
  
  //Rcpp::NumericMatrix Powers_dif(o1, m);
  int e2[o1];
  
  for ( j = 0; j < o1; j++ )
  {
    //f1 = mono_unrank_grlex ( m, e1[j] );
    Rcpp::IntegerVector f1 = Powers(j,Rcpp::_);
    for ( i = 0; i < m; i++ )
    {
      c2[j] = c2[j] * i4_fall ( f1[i], dif[i] );
      f1[i] = i4_max ( f1[i] - dif[i], 0 );
    }
    e2[j] = burkardt_rank_grlex ( m, f1);
    //delete [] f1;
  }
  
  polynomial_sort ( o2, c2, e2 );
  
  Rcpp::List out = burkardt_compress ( o2, c2, e2, m );
  return out;
}
//****************************************************************************80

double* burkardt_mono_value ( int m, int n, Rcpp::IntegerVector f, Rcpp::NumericMatrix x )
{
  int i;
  int j;
  double *v;
  
  v = new double[n];
  
  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
    for ( i = 0; i < m; i++ )
    {
      v[j] = v[j] * pow ( x[i+j*m], f[i] );
    }
  }
  
  return v;
}
//****************************************************************************80

// [[Rcpp::export]]
Rcpp::NumericVector burkardt_polynomial_value ( 
    Rcpp::NumericVector c, Rcpp::IntegerMatrix Powers,  
                           Rcpp::NumericMatrix x )
{
  int j;
  int k;
  Rcpp::NumericMatrix tx = Rcpp::transpose(x);
  int n = x.nrow();
  Rcpp::NumericVector p(n);
  int m = x.ncol();
  int o = c.size();
  
  for ( j = 0; j < o; j++ )
  {
    Rcpp::IntegerVector f = Powers(j, Rcpp::_);
    double* v = burkardt_mono_value ( m, n, f, tx );
    for ( k = 0; k < n; k++ )
    {
      p[k] = p[k] + c[j] * v[k];
    }
    delete[] v;
  }
  
  return p;
}


