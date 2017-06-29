#ifndef MY_VECTOR_H
#define MY_VECTOR_H

template <class T> class V
{
public:
	V<T> (T x_=0, T y_=0, T z_=0) {
		x=x_;
		y=y_;
		z=z_;
	}
	V<T> &operator =(V<T> v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	T x;
	T y;
	T z;
};


template <class T> V<T> operator +(V<T> v1, V<T> v2) {
	V<T> tmp(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
	return tmp;
}

template <class T> V<T> operator -(V<T> v1, V<T> v2) {
	V<T> tmp(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
	return tmp;
}

template <class T> T operator *(V<T> v1, V<T> v2) {
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

template <class T> V<T> operator *(double num, V<T> v) {
    V<T> tmp(v.x*num, v.y*num, v.z*num);
	return tmp;
}



#endif
