#pragma once
#include "includes.h"

template <class T>
class vector3
{
public:
	union{
		struct{T x,y,z;};
		struct{T coords[3];}; //ensure compatability with OpenGL
	};

	vector3(T _x, T _y, T _z)
	{
		x=_x;
		y=_y;
		z=_z;
	}
	
	vector3(T v[3])
	{
		x=v[0];
		y=v[1];
		z=v[2];
	}
	vector3(const vector3<T>& v)
	{
		x=v.x;
		y=v.y;
		z=v.z;
	}

	vector3(void){x=y=z=0;}


	inline T    dot(vector3<T>& v);
	inline vector3<T> cross(vector3<T>& v);
	inline vector3<T> unit();
	inline T    mag();

	inline vector3<T>& operator=(vector3<T>&); //ensure well defined assignment
	inline vector3<T> operator*(T);		//scalar mul
	inline vector3<T> operator*(vector3<T>& v);//scalar mul
	inline vector3<T> operator/(T);		//scalar divide
	inline vector3<T> operator+(vector3<T>& v);//add
	inline vector3<T>& operator+=(vector3<T>& v);//add
	inline vector3<T> operator-(vector3<T>& v);//subtract
	inline vector3<T>& operator-=(vector3<T>& v);//subtract
	inline vector3<T> operator-();			 //inverse
	inline vector3<T> operator~();			 //unit

	static vector3<T> lerp(vector3<T>& start, vector3<T>& end, T t);
	static vector3<T> nlerp(vector3<T>& start, vector3<T>& end, T t);
	static vector3<T> slerp(vector3<T>& start, vector3<T>& end, T t);
};

template <class T>
inline
vector3<T> vector3<T>::lerp(vector3<T>& start, vector3<T>& end, T t){
	return (start*(1-t) + end*t);
}

//assumes start and end are unit vectors
template <class T>
inline
vector3<T> vector3<T>::nlerp(vector3<T>& start, vector3<T>& end, T t){
	return lerp(start,end,t).unit();
}

//assumes start and end are unit vectors
template <class T>
inline
vector3<T> vector3<T>::slerp(vector3<T>& start, vector3<T>& end, T t){
	T innerp = start.dot(end);
	innerp=CLAMP(innerp,-1,1);
	T angle= acos(innerp) * t;
	vector3<T> base = (end - start * innerp).unit();
	return start * cos(angle) + base * sin(angle);
}

template <class T>
inline
vector3<T>& vector3<T>::operator=(vector3<T>& v){
	x=v.x;
	y=v.y;
	z=v.z;
	return *this;
}

template <class T>
inline
T vector3<T>::dot(vector3<T>& v){
	return x * v.x + y * v.y + z * v.z;
}

template <class T>
inline
vector3<T> vector3<T>::cross(vector3<T>& v){
	return vector3<T>(y * v.z - v.y * z, z * v.x - v.z * x, x * v.y - v.x * y);
}

template <class T>
inline
vector3<T> vector3<T>::unit(){
T m=mag();
return vector3<T>(x/m, y/m, z/m);
}

template <class T>
inline
vector3<T> vector3<T>::operator~(){
	return unit();
}

template <class T>
inline
vector3<T> vector3<T>::operator*(T k){
	return vector3<T>(x*k,y*k,z*k);
}

template <class T>
inline
vector3<T> operator*(T k, vector3<T>& v){
	return vector3<T>(v.x*k,v.y*k,v.z*k);
}

template <class T>
inline
vector3<T> vector3<T>::operator*(vector3<T>& v){
	return vector3<T>(x*v.x,y*v.y,z*v.z);
}

template <class T>
inline
vector3<T> vector3<T>::operator/(T k){
	return vector3<T>(x/k,y/k,z/k);
}

template <class T>
inline
vector3<T> vector3<T>::operator+(vector3<T>& v){
	return vector3<T>(x+v.x,y+v.y,z+v.z);
}

template <class T>
inline
vector3<T>& vector3<T>::operator+=(vector3<T>& v){
	*this = *this + v;
	return *this;
}

template <class T>
inline
vector3<T> vector3<T>::operator-(vector3<T>& v){
	return vector3<T>(x-v.x,y-v.y,z-v.z);
}

template <class T>
inline
vector3<T>& vector3<T>::operator-=(vector3<T>& v){
	*this = *this - v;
	return *this;
}

template <class T>
inline
vector3<T> vector3<T>::operator-(){
	return vector3<T>(-x, -y, -z);
}

template <class T>
inline
T vector3<T>::mag(){
	return sqrt(x*x+y*y+z*z);
}