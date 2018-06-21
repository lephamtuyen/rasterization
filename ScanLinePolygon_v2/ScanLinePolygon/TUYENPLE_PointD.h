#pragma once
#include <vector>

class PointD;

typedef std::vector<PointD> PointDArray;

class PointD {

public:
	double m_X;
	double m_Y;

	PointD(double x = 0.0, double y = 0.0);
	~PointD(void);

	double GetX() {
		return m_X;
	}

	double GetY() {
		return m_Y;
	}

	void SetXY(double x, double y) {
		m_X = x;
		m_Y = y;
	}

	void SetX(double x) {
		m_X = x;
	}

	void SetY(double y) {
		m_Y = y;
	}

	PointD& operator=( const PointD& rhs ) ;

	bool operator==(const PointD& rhs) const ;
	bool operator!=(const PointD& rhs) const ;

	PointD operator+(const PointD& rhs) const ;
	PointD& operator+=(const PointD& rhs) ;
	PointD operator-(const PointD& rhs) const ;
	PointD& operator-=(const PointD& rhs) ;
	PointD operator/( const double& dDenominator ) ;

	double GetDistanceWith(PointD Pt)  const;
};

