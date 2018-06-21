#pragma once

class RectD {

private:
	double m_left;
	double m_top;
	double m_right;
	double m_bottom;

public:
	RectD(double x1 = 0., double y1 = 0., double x2 = 0., double y2 = 0.);
	~RectD(void);

	RectD& operator=(const RectD& rhs) ;

	double GetLeft() {
		return m_left;
	}

	void SetLeft(double x) {
		m_left = x;	
	}

	double GetTop() {
		return m_top;
	}

	void SetTop(double y) {
		m_top = y;
	}

	double GetRight() {
		return m_right;
	}

	void SetRight(double x) {
		m_right = x;	
	}

	double GetBottom() {
		return m_bottom;
	}

	void SetBottom(double maxY) {
		m_bottom = maxY;
	}

};