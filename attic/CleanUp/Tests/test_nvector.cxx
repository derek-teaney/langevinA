#include "nvector.h"
#include <iostream>

int main(int argc, char **argv) {
  std::cout << "Hello World!" << std::endl;

  nvector<double,1> v1(4) ; 
  nvector<double,1> v2(8) ; 
  nvector<double,1> v3({12}) ; 
  nvector<double,1> v4({1}) ;
  v1.v = {1,2, 3, 4.2} ;
  std::cout << "Tests of v1 " << std::endl; 
  std::cout << v1.v.size()  << std::endl; 
  nvector<double,1> v5 ; 
  v1.print() ;
  for (int i=0; i < v1.N[0] ; i++) {
    std::cout << "v1 = " << v1(i)  << std::endl; 
  }
  
  v2.print() ;
  v3.print() ;
  v4.print() ;
  std::cout << "Tests of v1" << std::endl;
  std:: cout << v4(0) << std::endl; 
  v4.v = {3.14159} ; 
  std:: cout << v4(0) << std::endl; 
  v5.print() ;

  std::array<int, 2> NN = {5, 2};
  nvector<double, 2> a(NN[0], NN[1]);
  a.print();
  for (int i = 0; i < a.N[0]; i++) {
    for (int j = 0; j < a.N[1]; j++) {
      a(i, j) = i * 100 + j;
    }
  }


  for (int i = 0; i < a.N[0]; i++) {
    for (int j = 0; j < a.N[1]; j++) {
      std::cout << a(i, j) << " ";
    }
    std::cout << std::endl;
  }

  // Test pointer arithmetic
  double *aptr = &a(1,0) ;
  for (int j = 0 ; j < a.N[1] ; j++) {
     std::cout << a(1, j) << " compare to " << aptr[j] << ", ";
  }
  std::cout << std::endl;

  nvector<double, 1> b(5);
  b.print();
  b.resize(6);
  b.print();

  nvector<double, 2> c;
  c.print();
  c.resize(4, 2);
  c.print();

  nvector<double, 3> d(4, 5, 2);
  d.print();
  d.resize(4, 3, 2);
  for (int i = 0; i < d.N[0]; i++) {
    for (int j = 0; j < d.N[1]; j++) {
      for (int k = 0; k < d.N[2]; k++) {
        d(i, j, k) = 10000 * i + 100 * j + k;
      }
    }
  }
  for (int i = 0; i < d.N[0]; i++) {
    for (int j = 0; j < d.N[1]; j++) {
      for (int k = 0; k < d.N[2]; k++) {
        std::cout << d(i, j, k) << " ";
      }
    }
    std::cout << std::endl;
  }

  int i=4 ;
  int j=2 ;
  nvector<double, 4> e(i, 5, 3, j);
  e.print();
  e.resize(i, 3, j, 1);
  e.print();

  nvector<double,5> f(4, 5, 3, 3, 2);
  f.print();
  f.resize(4, 6, 2, 3, 2);
  f.print();
}
