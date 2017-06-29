#include <iostream>
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"

boost::numeric::ublas::matrix<double> operator *(boost::numeric::ublas::matrix<double> m1, boost::numeric::ublas::matrix<double> m2) {
	if (m1.size2() != m2.size1())
		cerr << "error: in the matrix<double> operator*(matrix<double> m1, matrix<double> m2);  matrix sizes do not match each other" << endl;
	boost::numeric::ublas::matrix<double> tmp(m1.size1(),m2.size2());
	for (unsigned i=0; i<tmp.size1(); ++i)
		for (unsigned j=0; j<tmp.size2(); ++j) 
			tmp(i,j)=0;
	for (unsigned i1=0; i1<m1.size1(); ++i1)
        for (unsigned j2=0; j2<m2.size2(); ++j2)
			for (unsigned j1=0; j1<m1.size2(); ++j1)
				tmp(i1,j2) += m1(i1,j1)*m2(j1,j2);
	return tmp;
}


boost::numeric::ublas::vector<double> operator *(boost::numeric::ublas::vector<double> v, boost::numeric::ublas::matrix<double> m) {
	if (v.size() != m.size1())
		cerr << "error: in the vector<double> operator *(vector<double> v, matrix<double> m);  matrix and vector sizes do not match each other" << endl;
	boost::numeric::ublas::vector<double> tmp(m.size2());
	for (unsigned i=0; i<tmp.size(); ++i) tmp(i)=0;
	cout << tmp << endl;
	for (unsigned i=0; i<m.size2(); ++i)
		for (unsigned j=0; j<v.size(); ++j)
			tmp(i) += v(j)*m(j,i);
	return tmp;
}

boost::numeric::ublas::vector<double> operator *(boost::numeric::ublas::matrix<double> m, boost::numeric::ublas::vector<double> v) {
	if (m.size2() != v.size())
		cerr << "error: in the vector<double> operator *(matrix<double> m, vector<double> v);  matrix and vector sizes do not match each other" << endl;
	boost::numeric::ublas::vector<double> tmp(m.size1());
	for (unsigned i=0; i<tmp.size(); ++i) tmp(i)=0;
	cout << tmp << endl;
	for (unsigned i=0; i<m.size1(); ++i)
		for (unsigned j=0; j<v.size(); ++j)
			tmp(i) += v(j)*m(i,j);
	return tmp;
}

boost::numeric::ublas::matrix<double> Transpose(boost::numeric::ublas::matrix<double> m) {
	boost::numeric::ublas::matrix<double> tmp(m.size1(),m.size2());
	for(int i=0; i<m.size1(); i++)
		for(int j=0; j<m.size2(); j++)
			tmp(i,j) = m(j,i);
	return tmp;
}

