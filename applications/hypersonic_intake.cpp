#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
using namespace std;

//--------------------------------------------------------------------------------
//>GLOBAL VARIABLES
//--------------------------------------------------------------------------------
const int cnum = 33600;
const int nnum = 38194;
int connect_num = 192712;

const int snum = 6;
const int vnum = 8;
const int rnum = 1 + snum;
#define PI acos(-1)
const int pnum = 5;
const int mnnum = 7;
const int mtnum = 5;
const int minum = 3;

double cfl;
double sim_time;
double max_time;
double eps;
double res;
double dt;
int iter;
int output;

int ck;
int ct;
double Pr;
double gam;
double mu_ref;

double dens[rnum * cnum];
double xmom[rnum * cnum];
double ymom[rnum * cnum];
double zmom[rnum * cnum];
double energy[rnum * cnum];

int element[vnum * cnum];
int connect[snum * cnum];
double notes[3 * nnum];
double center[3 * cnum];
double normal[3 * snum * cnum];
double points[3 * snum * cnum];
double scale[(1 + snum) * cnum];

double rho_int, u_int, v_int, w_int, lam_int, p_int;
double rho_in, u_in, v_in, w_in, lam_in, p_in;
double rho_out, u_out, v_out, w_out, lam_out, p_out;

//--------------------------------------------------------------------------------
//>COMMENLY USED FUNCTION
//--------------------------------------------------------------------------------

//>function: obtain ratio of specific heat
double get_gamma() {
	return (double)(ck + ct + 2) / (ck + ct);
}

//>function: obtain primary variables
void get_primary(double prim[pnum], double w[pnum]) {
	prim[0] = w[0];
	prim[1] = w[1] / w[0];
	prim[2] = w[2] / w[0];
	prim[3] = w[3] / w[0];
	prim[4] = 0.5 * w[0] / (gam - 1.0) / (w[4] - 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]) / w[0]);
}

//>function: obtain conservative variables
void get_conserved(double w[pnum], double prim[pnum]) {
	w[0] = prim[0];
	w[1] = prim[0] * prim[1];
	w[2] = prim[0] * prim[2];
	w[3] = prim[0] * prim[3];
	w[4] = 0.5 * prim[0] / prim[4] / (gam - 1.0) + 0.5 * prim[0] * (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3]);
}

//>function: frame rotation
void rotation(double w[pnum], double n1, double n2, double n3, int drc) {

	double a1, a2, a3, b1, b2, b3, c1, c2, c3;
	double o0 = w[0], o1 = w[1], o2 = w[2], o3 = w[3], o4 = w[4];

	if (n1 == -1) {
		a1 = -1.0; a2 = 0.0;  a3 = 0.0;
		b1 = 0.0;  b2 = -1.0; b3 = 0.0;
		c1 = 0.0;  c2 = 0.0;  c3 = 1.0;
	}
	else {
		a1 = n1;  a2 = n2;                a3 = n3;
		b1 = -n2; b2 = n1 + n3 * n3 / (1.0 + n1); b3 = -n2 * n3 / (1.0 + n1);
		c1 = -n3; c2 = -n2 * n3 / (1.0 + n1);   c3 = 1.0 - n3 * n3 / (1.0 + n1);
	}

	if (drc == -1) {
		double ar, ar1, ar2, ar3, br1, br2, br3, cr1, cr2, cr3;
		ar1 = a1; br1 = a2; cr1 = a3;
		ar2 = b1; br2 = b2; cr2 = b3;
		ar3 = c1; br3 = c2; cr3 = c3;

		ar = ar1 * (br2 * cr3 - cr2 * br3) - ar2 * (br1 * cr3 - cr1 * br3) + ar3 * (br1 * cr2 - cr1 * br2);

		a1 = (br2 * cr3 - cr2 * br3) / ar; a2 = (cr1 * br3 - br1 * cr3) / ar; a3 = (br1 * cr2 - cr1 * br2) / ar;
		b1 = (cr2 * ar3 - ar2 * cr3) / ar; b2 = (ar1 * cr3 - cr1 * ar3) / ar; b3 = (ar2 * cr1 - ar1 * cr2) / ar;
		c1 = (ar2 * br3 - br2 * ar3) / ar; c2 = (br1 * ar3 - ar1 * br3) / ar; c3 = (ar1 * br2 - ar2 * br1) / ar;
	}

	w[0] = o0;
	w[1] = a1 * o1 + a2 * o2 + a3 * o3;
	w[2] = b1 * o1 + b2 * o2 + b3 * o3;
	w[3] = c1 * o1 + c2 * o2 + c3 * o3;
	w[4] = o4;
}

//>function: SIGN function
int sign(double x) {
	if (x > 0) { return 1; }
	else if (x == 0) { return 0; }
	else { return -1; }
}

//>function: print 1D array
template <typename type1>
void print1(type1 array, int size) {
	for (int i = 0; i < size; i++) {
		cout << array[i] << ", ";
	}
	cout << endl;
}

//--------------------------------------------------------------------------------
//>MESH INPUT
//--------------------------------------------------------------------------------

//>function: read mesh file
void read_mesh() {
	ifstream infile;
	infile.open("mesh/hypersonic_intake/notes.txt");
	for (int i = 0; i < 3 * nnum; i++) {
		infile >> notes[i];
	}
	infile.close();

	infile.open("mesh/hypersonic_intake/elements.txt");
	for (int i = 0; i < vnum * cnum; i++) {
		infile >> element[i];
		element[i] = element[i] - 1;
	}
	infile.close();

	for (int i = 0; i < snum * cnum; i++) { connect[i] = -1; }
	infile.open("mesh/hypersonic_intake/connectivity.txt");
	int cell, side, near;
	for (int i = 0; i < connect_num; i++) {
		infile >> cell;
		infile >> side;
		infile >> near;
		connect[snum * (cell - 1) + side - 1] = near - 1;
	}
	infile.close();
}

//>function: non-zero zone
void no_zero() {
	for (int i = 0; i < 3 * nnum; i++) {
		if (notes[i] < 0.00005) { notes[i] = 0; }
	}
}

//>function: geometry information of face 
double face(double n[3], double c[3], double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4) {

	double a1, a2, a3, b1, b2, b3, n1, n2, n3, m;

	a1 = x3 - x1; a2 = y3 - y1; a3 = z3 - z1;
	b1 = x4 - x1; b2 = y4 - y1; b3 = z4 - z1;

	n1 = a2 * b3 - a3 * b2;
	n2 = a3 * b1 - a1 * b3;
	n3 = a1 * b2 - a2 * b1;

	m = sqrt(n1 * n1 + n2 * n2 + n3 * n3);

	n[0] = n1 / m;
	n[1] = n2 / m;
	n[2] = n3 / m;

	double t1, t2, t3, t;

	t1 = ((y3 - y1) / (y2 - y1) - (x3 - x1) / (x2 - x1)) / ((x4 - x3) / (x2 - x1) - (y4 - y3) / (y2 - y1));
	t2 = ((z3 - z1) / (z2 - z1) - (x3 - x1) / (x2 - x1)) / ((x4 - x3) / (x2 - x1) - (z4 - z3) / (z2 - z1));
	t3 = ((z3 - z1) / (z2 - z1) - (y3 - y1) / (y2 - y1)) / ((y4 - y3) / (y2 - y1) - (z4 - z3) / (z2 - z1));

	if (!isnan(t1)) { t = t1; }
	else if (!isnan(t2)) { t = t2; }
	else { t = t3; }

	c[0] = x3 + (x4 - x3) * t;
	c[1] = y3 + (y4 - y3) * t;
	c[2] = z3 + (z4 - z3) * t;

	double d1, d2, d3, d4, d5, k1, k2;

	d1 = sqrt(pow(x3 - x1, 2) + pow(y3 - y1, 2) + pow(z3 - z1, 2));
	d2 = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2) + pow(z3 - z2, 2));
	d3 = sqrt(pow(x4 - x2, 2) + pow(y4 - y2, 2) + pow(z4 - z2, 2));
	d4 = sqrt(pow(x1 - x4, 2) + pow(y1 - y4, 2) + pow(z1 - z4, 2));
	d5 = sqrt(pow(x3 - x4, 2) + pow(y3 - y4, 2) + pow(z3 - z4, 2));
	k1 = 0.5 * (d1 + d4 + d5);
	k2 = 0.5 * (d2 + d3 + d5);
	return sqrt(k1 * (k1 - d1) * (k1 - d4) * (k1 - d5)) + sqrt(k2 * (k2 - d2) * (k2 - d3) * (k2 - d5));

}

//>function: geometry information (3D Hex)
void geometry_hex() {

	int n[vnum], nn[vnum];
	double x[vnum], y[vnum], z[vnum];
	int n1, n2, n3, n4, nr, nrr;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	double a, crs[3], nor[3], cen[3];
	double cz1, cz2;
	double s1, s2, s3, s4, s5, s6;

	for (int i = 0; i < cnum; i++) {
		for (int j = 0; j < vnum; j++) {
			n[j] = element[i * vnum + j];
			x[j] = notes[3 * n[j]];
			y[j] = notes[3 * n[j] + 1];
			z[j] = notes[3 * n[j] + 2];
		}
		n1 = 0; n2 = 2; n3 = 1; n4 = 3; nr = 0; nrr = 0;
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);  if (nor[1] != 0) { cout << i << endl; }
		for (int j = 0; j < 4; j++) {
			if (x[j] > crs[0] && y[j] > crs[1]) { nr = nr + 1; if (nr > 1) { if (y[nn[0]] < y[j]) { nn[3] = nn[0]; nn[0] = j; } else { nn[3] = j; } } else { nn[0] = j; } }
			else if (x[j] < crs[0] && y[j] < crs[1]) { nrr = nrr + 1; if (nrr > 1) { if (y[nn[2]] > y[j]) { nn[1] = nn[2]; nn[2] = j; } else { nn[1] = j; } } else { nn[2] = j; } }
			else if (x[j] > crs[0]) { nn[3] = j; }
			else { nn[1] = j; }
		}

		n1 = nn[0]; n2 = nn[2]; n3 = nn[1]; n4 = nn[3];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		cz1 = crs[2];

		n1 = nn[0] + 4; n2 = nn[2] + 4; n3 = nn[1] + 4; n4 = nn[3] + 4;
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		cz2 = crs[2];

		for (int j = 0; j < 4; j++) {
			if (cz1 > cz2) { nn[j + 4] = nn[j] + 4; }
			else { nn[j + 4] = nn[j]; nn[j] = nn[j] + 4; }
		}
		for (int j = 0; j < vnum; j++) { element[i * vnum + j] = n[nn[j]]; }


		n1 = nn[0]; n2 = nn[6]; n3 = nn[2]; n4 = nn[4];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, cen, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		center[3 * i] = cen[0];
		center[3 * i + 1] = cen[1];
		center[3 * i + 2] = cen[2];

		n1 = nn[2]; n2 = nn[5]; n3 = nn[1]; n4 = nn[6];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 1] = a;
		normal[3 * snum * i] = nor[0];
		normal[3 * snum * i + 1] = nor[1];
		normal[3 * snum * i + 2] = nor[2];
		points[3 * snum * i] = crs[0];
		points[3 * snum * i + 1] = crs[1];
		points[3 * snum * i + 2] = crs[2];
		s1 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		n1 = nn[0]; n2 = nn[7]; n3 = nn[3]; n4 = nn[4];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 2] = a;
		normal[3 * snum * i + 3] = nor[0];
		normal[3 * snum * i + 4] = nor[1];
		normal[3 * snum * i + 5] = nor[2];
		points[3 * snum * i + 3] = crs[0];
		points[3 * snum * i + 4] = crs[1];
		points[3 * snum * i + 5] = crs[2];
		s2 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		n1 = nn[3]; n2 = nn[6]; n3 = nn[2]; n4 = nn[7];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 3] = a;
		normal[3 * snum * i + 6] = nor[0];
		normal[3 * snum * i + 7] = nor[1];
		normal[3 * snum * i + 8] = nor[2];
		points[3 * snum * i + 6] = crs[0];
		points[3 * snum * i + 7] = crs[1];
		points[3 * snum * i + 8] = crs[2];
		s3 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		n1 = nn[1]; n2 = nn[4]; n3 = nn[0]; n4 = nn[5];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 4] = a;
		normal[3 * snum * i + 9] = nor[0];
		normal[3 * snum * i + 10] = nor[1];
		normal[3 * snum * i + 11] = nor[2];
		points[3 * snum * i + 9] = crs[0];
		points[3 * snum * i + 10] = crs[1];
		points[3 * snum * i + 11] = crs[2];
		s4 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		n1 = nn[7]; n2 = nn[5]; n3 = nn[6]; n4 = nn[4];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 5] = a;
		normal[3 * snum * i + 12] = nor[0];
		normal[3 * snum * i + 13] = nor[1];
		normal[3 * snum * i + 14] = nor[2];
		points[3 * snum * i + 12] = crs[0];
		points[3 * snum * i + 13] = crs[1];
		points[3 * snum * i + 14] = crs[2];
		s5 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		n1 = nn[0]; n2 = nn[2]; n3 = nn[1]; n4 = nn[3];
		x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4];
		y1 = y[n1]; y2 = y[n2]; y3 = y[n3]; y4 = y[n4];
		z1 = z[n1]; z2 = z[n2]; z3 = z[n3]; z4 = z[n4];
		a = face(nor, crs, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
		scale[rnum * i + 6] = a;
		normal[3 * snum * i + 15] = nor[0];
		normal[3 * snum * i + 16] = nor[1];
		normal[3 * snum * i + 17] = nor[2];
		points[3 * snum * i + 15] = crs[0];
		points[3 * snum * i + 16] = crs[1];
		points[3 * snum * i + 17] = crs[2];
		s6 = a / 3 * fabs(nor[0] * cen[0] + nor[1] * cen[1] + nor[2] * cen[2] - nor[0] * x1 - nor[1] * y1 - nor[2] * z1);

		scale[rnum * i] = s1 + s2 + s3 + s4 + s5 + s6;
	}

	int k, cnew[snum];
	double xj, yj, zj, xk, yk, zk;
	double smin, s[snum];

	for (int i = 0; i < cnum; i++) {
		for (int j = 0; j < snum; j++) { cnew[j] = -1; }
		for (int j = 0; j < snum; j++) {

			k = connect[snum * i + j];
			smin = 10000;
			if (k != -1) {
				xk = center[3 * k]; yk = center[3 * k + 1]; zk = center[3 * k + 2];

				for (int m = 0; m < snum; m++) {
					xj = points[3 * snum * i + 3 * m]; yj = points[3 * snum * i + 3 * m + 1]; zj = points[3 * snum * i + 3 * m + 2];
					s[m] = sqrt(pow(xj - xk, 2) + pow(yj - yk, 2) + pow(zj - zk, 2));
					smin = min(smin, s[m]);
				}
				for (int m = 0; m < snum; m++) {
					if (smin == s[m]) { cnew[m] = k; break; }
				}
			}

		}
		for (int j = 0; j < snum; j++) { connect[snum * i + j] = cnew[j]; }
	}

}

//--------------------------------------------------------------------------------
//>FIND BOUNDARY
//--------------------------------------------------------------------------------

//>function: find boundary sides
int find_bc(int bc[], int size, double value, int drc) {

	for (int i = 0; i < size; i++) { bc[i] = -1; }
	int num = 0;

	for (int i = 0; i < cnum; i++) {
		for (int j = 0; j < snum; j++) {
			if (drc == 0) {
				if (points[3 * snum * i + 3 * j] == value) { bc[2 * num] = i; bc[2 * num + 1] = j; num++; }
			}
			else if (drc == 1) {
				if (points[3 * snum * i + 3 * j + 1] == value) { bc[2 * num] = i; bc[2 * num + 1] = j; num++; }
			}
			else {
				if (points[3 * snum * i + 3 * j + 2] == value) { bc[2 * num] = i; bc[2 * num + 1] = j; num++; }
			}
		}
	}

	return num;
}

//>function: implementing periodic connection
void connect_periodic() {

	const int size = 1000;
	int bc0[size], bc1[size], bc2[size], bc3[size], bc4[size], bc5[size];
	int n0, n1, n2, n3, n4, n5;
	int cell, side, cellr, sider;
	double x, y, z, xr, yr, zr;

	n0 = find_bc(bc0, size, 0, 0);
	n1 = find_bc(bc1, size, 2, 0);
	n2 = find_bc(bc2, size, 0, 1);
	n3 = find_bc(bc3, size, 2, 1);
	n4 = find_bc(bc4, size, 0, 2);
	n5 = find_bc(bc5, size, 2, 2);

	for (int i = 0; i < n0; i++) {
		cell = bc0[2 * i];
		side = bc0[2 * i + 1];
		y = points[3 * snum * cell + 3 * side + 1];
		z = points[3 * snum * cell + 3 * side + 2];

		for (int j = 0; j < n1; j++) {
			cellr = bc1[2 * j];
			sider = bc1[2 * j + 1];
			yr = points[3 * snum * cellr + 3 * sider + 1];
			zr = points[3 * snum * cellr + 3 * sider + 2];

			if (y == yr && z == zr) {
				connect[snum * cell + side] = cellr;
				connect[snum * cellr + sider] = cell;
			}
		}
	}

	for (int i = 0; i < n2; i++) {
		cell = bc2[2 * i];
		side = bc2[2 * i + 1];
		x = points[3 * snum * cell + 3 * side];
		z = points[3 * snum * cell + 3 * side + 2];

		for (int j = 0; j < n3; j++) {
			cellr = bc3[2 * j];
			sider = bc3[2 * j + 1];
			xr = points[3 * snum * cellr + 3 * sider];
			zr = points[3 * snum * cellr + 3 * sider + 2];

			if (x == xr && z == zr) {
				connect[snum * cell + side] = cellr;
				connect[snum * cellr + sider] = cell;
			}
		}
	}

	for (int i = 0; i < n4; i++) {
		cell = bc4[2 * i];
		side = bc4[2 * i + 1];
		y = points[3 * snum * cell + 3 * side + 1];
		x = points[3 * snum * cell + 3 * side];

		for (int j = 0; j < n5; j++) {
			cellr = bc5[2 * j];
			sider = bc5[2 * j + 1];
			yr = points[3 * snum * cellr + 3 * sider + 1];
			xr = points[3 * snum * cellr + 3 * sider];

			if (y == yr && x == xr) {
				connect[snum * cell + side] = cellr;
				connect[snum * cellr + sider] = cell;
			}
		}
	}

}

//--------------------------------------------------------------------------------
//>RECONSTRUCTION
//--------------------------------------------------------------------------------

//>function: First-order reconstruction
void reconstruct_1st() {
	for (int i = 0; i < cnum; i++) {
		for (int j = 1; j < rnum; j++) {
			dens[rnum * i + j] = dens[rnum * i];
			xmom[rnum * i + j] = xmom[rnum * i];
			ymom[rnum * i + j] = ymom[rnum * i];
			zmom[rnum * i + j] = zmom[rnum * i];
			energy[rnum * i + j] = energy[rnum * i];
		}
	}
}

//>function: Second-order reconstruction for quad cell
void reconstruct_2nd_quad() {

	int near0, near1, near2, near3;
	double d0, d1, d2, d3, k0, k1, k2, k3;
	double wc[pnum], w0[pnum], w1[pnum], w2[pnum], w3[pnum];
	double s0, s1, s2, s3, sx, sy;

	for (int i = 0; i < cnum; i++) {

		near0 = connect[snum * i];    if (near0 == -1) { near0 = i; }
		near1 = connect[snum * i + 1];  if (near1 == -1) { near1 = i; }
		near2 = connect[snum * i + 2];  if (near2 == -1) { near2 = i; }
		near3 = connect[snum * i + 3];  if (near3 == -1) { near3 = i; }

		d0 = points[3 * snum * i + 1] - center[3 * i + 1];   k0 = center[3 * near0 + 1] - center[3 * i + 1];   if (fabs(k0) > 5 * d0) { k0 = 2 * d0; } //-------only for this periodic condition
		d1 = center[3 * i] - points[3 * snum * i + 3];     k1 = center[3 * i] - center[3 * near1];       if (fabs(k1) > 5 * d1) { k1 = 2 * d1; }
		d2 = center[3 * i + 1] - points[3 * snum * i + 7];   k2 = center[3 * i + 1] - center[3 * near2 + 1];   if (fabs(k2) > 5 * d2) { k2 = 2 * d2; }
		d3 = points[3 * snum * i + 9] - center[3 * i];     k3 = center[3 * near3] - center[3 * i];       if (fabs(k3) > 5 * d3) { k3 = 2 * d3; }

		wc[0] = dens[rnum * i];
		wc[1] = xmom[rnum * i];
		wc[2] = ymom[rnum * i];
		wc[3] = zmom[rnum * i];
		wc[4] = energy[rnum * i];

		w0[0] = dens[rnum * near0];
		w0[1] = xmom[rnum * near0];
		w0[2] = ymom[rnum * near0];
		w0[3] = zmom[rnum * near0];
		w0[4] = energy[rnum * near0];

		w1[0] = dens[rnum * near1];
		w1[1] = xmom[rnum * near1];
		w1[2] = ymom[rnum * near1];
		w1[3] = zmom[rnum * near1];
		w1[4] = energy[rnum * near1];

		w2[0] = dens[rnum * near2];
		w2[1] = xmom[rnum * near2];
		w2[2] = ymom[rnum * near2];
		w2[3] = zmom[rnum * near2];
		w2[4] = energy[rnum * near2];

		w3[0] = dens[rnum * near3];
		w3[1] = xmom[rnum * near3];
		w3[2] = ymom[rnum * near3];
		w3[3] = zmom[rnum * near3];
		w3[4] = energy[rnum * near3];

		for (int j = 0; j < pnum; j++) {
			s0 = (w0[j] - wc[j]) / k0;
			s2 = (wc[j] - w2[j]) / k2;
			if (near0 == i) { sy = s2; }
			else if (near2 == i) { sy = s0; }
			else { sy = (sign(s0) + sign(s2)) * fabs(s0 * s2) / (fabs(s0) + fabs(s2)); if (isnan(sy)) { sy = 0; } }

			s1 = (wc[j] - w1[j]) / k1;
			s3 = (w3[j] - wc[j]) / k3;
			if (near1 == i) { sx = s3; }
			else if (near3 == i) { sx = s1; }
			else { sx = (sign(s1) + sign(s3)) * fabs(s1 * s3) / (fabs(s1) + fabs(s3)); if (isnan(sx)) { sx = 0; } }

			if (j == 0) {
				dens[rnum * i + 1] = wc[j] + sy * d0;
				dens[rnum * i + 2] = wc[j] - sx * d1;
				dens[rnum * i + 3] = wc[j] - sy * d2;
				dens[rnum * i + 4] = wc[j] + sx * d3;
			}
			else if (j == 1) {
				xmom[rnum * i + 1] = wc[j] + sy * d0;
				xmom[rnum * i + 2] = wc[j] - sx * d1;
				xmom[rnum * i + 3] = wc[j] - sy * d2;
				xmom[rnum * i + 4] = wc[j] + sx * d3;
			}
			else if (j == 2) {
				ymom[rnum * i + 1] = wc[j] + sy * d0;
				ymom[rnum * i + 2] = wc[j] - sx * d1;
				ymom[rnum * i + 3] = wc[j] - sy * d2;
				ymom[rnum * i + 4] = wc[j] + sx * d3;
			}
			else if (j == 3) {
				zmom[rnum * i + 1] = wc[j] + sy * d0;
				zmom[rnum * i + 2] = wc[j] - sx * d1;
				zmom[rnum * i + 3] = wc[j] - sy * d2;
				zmom[rnum * i + 4] = wc[j] + sx * d3;
			}
			else {
				energy[rnum * i + 1] = wc[j] + sy * d0;
				energy[rnum * i + 2] = wc[j] - sx * d1;
				energy[rnum * i + 3] = wc[j] - sy * d2;
				energy[rnum * i + 4] = wc[j] + sx * d3;
			}
		}

	}

}

//--------------------------------------------------------------------------------
//>GKS FLUX CALCULATION
//--------------------------------------------------------------------------------

//>function: obtain collision time
double get_tao(double lamo, double rhoo, double pl, double pr) {
	double c1 = 0.05, c2 = 5.0, c3 = 0.0;
	double tao1, tao2, tao3;
	tao1 = c1 * dt;
	tao2 = c2 * fabs((pl - pr) / (pl + pr)) * dt;
	tao3 = c3 * mu_ref * 2 * lamo / rhoo;
	return tao1 + tao2 + tao3;
}
//inviscid: c1=0.05 c2=5.0 c3=0.0
//non-equilibrium: c1=0.0 c2=1.0 c3=1.0

//>function: obtain macroscopic slope
void get_slope(double a[pnum], double prim[pnum], double s[pnum], int sg) {
	s[0] = sg * s[0] / prim[0]; s[1] = sg * s[1] / prim[0]; s[2] = sg * s[2] / prim[0]; s[3] = sg * s[3] / prim[0]; s[4] = sg * s[4] / prim[0];

	double r4 = 2 * s[4] - (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3] + 0.5 * (ck + ct) / prim[4]) * s[0];
	double r3 = s[3] - prim[3] * s[0];
	double r2 = s[2] - prim[2] * s[0];
	double r1 = s[1] - prim[1] * s[0];

	a[4] = 4.0 * prim[4] * prim[4] / (ck + ct) * (r4 - 2 * prim[1] * r1 - 2 * prim[2] * r2 - 2 * prim[3] * r3);
	a[3] = 2 * prim[4] * r3 - prim[3] * a[4];
	a[2] = 2 * prim[4] * r2 - prim[2] * a[4];
	a[1] = 2 * prim[4] * r1 - prim[1] * a[4];
	a[0] = s[0] - prim[1] * a[1] - prim[2] * a[2] - prim[3] * a[3] - 0.5 * a[4] * (prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3] + 0.5 * (ck + ct) / prim[4]);
}

//>function£ºobtain velocity moments
void moment_u(double mu[mnnum], double mul[mnnum], double mur[mnnum], double mv[mtnum], double mw[mtnum], double mi[minum], double prim[pnum]) {
	mul[0] = 0.5 * erfc(-sqrt(prim[4]) * prim[1]);
	mul[1] = prim[1] * mul[0] + 0.5 * exp(-prim[4] * prim[1] * prim[1]) / sqrt(PI * prim[4]);
	mur[0] = 0.5 * erfc(sqrt(prim[4]) * prim[1]);
	mur[1] = prim[1] * mur[0] - 0.5 * exp(-prim[4] * prim[1] * prim[1]) / sqrt(PI * prim[4]);
	mu[0] = 1.0; mu[1] = prim[1];
	mv[0] = 1.0; mv[1] = prim[2];
	mw[0] = 1.0; mw[1] = prim[3];

	for (int i = 2; i < mnnum; i++) {
		mul[i] = prim[1] * mul[i - 1] + 0.5 * (i - 1) * mul[i - 2] / prim[4];
		mur[i] = prim[1] * mur[i - 1] + 0.5 * (i - 1) * mur[i - 2] / prim[4];
		mu[i] = mul[i] + mur[i];
		if (i < mtnum) {
			mv[i] = prim[2] * mv[i - 1] + 0.5 * (i - 1) * mv[i - 2] / prim[4];
			mw[i] = prim[3] * mw[i - 1] + 0.5 * (i - 1) * mw[i - 2] / prim[4];
		}
	}

	mi[0] = 1.0;
	mi[1] = 0.5 * ck / prim[4];
	mi[2] = (ck * ck + 2.0 * ck) / (4.0 * prim[4] * prim[4]);
}

void moment_uvw(double muvw[pnum], int ua, int ub, int uc, int ud, double mu[mnnum], double mv[mtnum], double mw[mtnum], double mi[minum], double prim[pnum]) {
	muvw[0] = prim[0] * mu[ua] * mv[ub] * mw[uc] * mi[ud / 2];
	muvw[1] = prim[0] * mu[ua + 1] * mv[ub] * mw[uc] * mi[ud / 2];
	muvw[2] = prim[0] * mu[ua] * mv[ub + 1] * mw[uc] * mi[ud / 2];
	muvw[3] = prim[0] * mu[ua] * mv[ub] * mw[uc + 1] * mi[ud / 2];
	muvw[4] = prim[0] * 0.5 * (mu[ua + 2] * mv[ub] * mw[uc] * mi[ud / 2] + mu[ua] * mv[ub + 2] * mw[uc] * mi[ud / 2] + mu[ua] * mv[ub] * mw[uc + 2] * mi[ud / 2] + mu[ua] * mv[ub] * mw[uc] * mi[(ud + 2) / 2]);
}

void moment_a(double ma[pnum], double a[pnum], int ua1, int ua2, int ua3, double mu[mnnum], double mv[mtnum], double mw[mtnum], double mi[minum], double prim[pnum]) {
	double muvw0[pnum], muvw1[pnum], muvw2[pnum], muvw3[pnum], muvw4[pnum], muvw5[pnum], muvw6[pnum], muvw7[pnum];
	moment_uvw(muvw0, ua1 + 0, ua2 + 0, ua3 + 0, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw1, ua1 + 1, ua2 + 0, ua3 + 0, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw2, ua1 + 0, ua2 + 1, ua3 + 0, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw3, ua1 + 0, ua2 + 0, ua3 + 1, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw4, ua1 + 2, ua2 + 0, ua3 + 0, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw5, ua1 + 0, ua2 + 2, ua3 + 0, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw6, ua1 + 0, ua2 + 0, ua3 + 2, 0, mu, mv, mw, mi, prim);
	moment_uvw(muvw7, ua1 + 0, ua2 + 0, ua3 + 0, 2, mu, mv, mw, mi, prim);

	for (int i = 0; i < pnum; i++) {
		ma[i] = a[0] * muvw0[i] + a[1] * muvw1[i] + a[2] * muvw2[i] + a[3] * muvw3[i] + 0.5 * a[4] * (muvw4[i] + muvw5[i] + muvw6[i] + muvw7[i]);
	}
}

double mompr(double a[pnum], int pr1, int pr2, int pr3, int pr4, double mu[mnnum], double mv[mtnum], double mw[mtnum], double mi[minum], double prim[pnum]) {
	double m0, m1, m2, m3, m4;
	m0 = prim[0] * mu[pr1] * mv[pr2] * mw[pr3] * mi[pr4 / 2];
	m1 = prim[0] * mu[pr1 + 1] * mv[pr2] * mw[pr3] * mi[pr4 / 2];
	m2 = prim[0] * mu[pr1] * mv[pr2 + 1] * mw[pr3] * mi[pr4 / 2];
	m3 = prim[0] * mu[pr1] * mv[pr2] * mw[pr3 + 1] * mi[pr4 / 2];
	m4 = prim[0] * 0.5 * (mu[pr1 + 2] * mv[pr2] * mw[pr3] * mi[pr4 / 2] + mu[pr1] * mv[pr2 + 2] * mw[pr3] * mi[pr4 / 2] + mu[pr1] * mv[pr2] * mw[pr3 + 2] * mi[pr4 / 2] + mu[pr1] * mv[pr2] * mw[pr3] * mi[(pr4 + 2) / 2]);
	return a[0] * m0 + a[1] * m1 + a[2] * m2 + a[3] * m3 + a[4] * m4;
}

//>main function
void fluxGKS(double flux[pnum], double wcl[pnum], double wcr[pnum], double wl[pnum], double wr[pnum], double dil, double dir) {

	double priml[pnum], primr[pnum], primo[pnum], wo[pnum];
	double mu[mnnum], mul[mnnum], mur[mnnum], mv[mtnum], mw[mtnum], mi[minum];
	double mu_l[mnnum], mul_l[mnnum], mur_l[mnnum], mv_l[mtnum], mw_l[mtnum], mi_l[minum];
	double mu_r[mnnum], mul_r[mnnum], mur_r[mnnum], mv_r[mtnum], mw_r[mtnum], mi_r[minum];
	double mom_o[pnum], mom_l[pnum], mom_r[pnum], mom_al[pnum], mom_ar[pnum], mom_aL[pnum], mom_aR[pnum], mom_Al[pnum], mom_Ar[pnum], mom_A[pnum];
	double w_l[pnum], w_r[pnum], w_L[pnum], w_R[pnum], w_T[pnum];
	double a_l[pnum], a_r[pnum], a_L[pnum], a_R[pnum], A_l[pnum], A_r[pnum], A[pnum];
	double tao, t0, t1, t2, t3, t4, t5;
	double lam_wo, rho_wo, p_wl, p_wr;

	//>equilibrium state at the interface
	get_primary(priml, wl);
	get_primary(primr, wr);

	moment_u(mu_l, mul_l, mur_l, mv_l, mw_l, mi_l, priml);
	moment_u(mu_r, mul_r, mur_r, mv_r, mw_r, mi_r, primr);
	moment_uvw(mom_l, 0, 0, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_uvw(mom_r, 0, 0, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);

	for (int i = 0; i < pnum; i++) { wo[i] = mom_l[i] + mom_r[i]; }
	get_primary(primo, wo);
	moment_u(mu, mul, mur, mv, mw, mi, primo);

	//>collision time and time integration terms
	lam_wo = primo[4]; rho_wo = primo[0]; p_wl = 0.5 * priml[0] / priml[4]; p_wr = 0.5 * primr[0] / primr[4];
	tao = get_tao(lam_wo, rho_wo, p_wl, p_wr);

	t0 = dt - tao * (1 - exp(-dt / tao));
	t1 = -(1 - exp(-dt / tao)) / t0;
	t2 = (-dt + 2 * tao * (1 - exp(-dt / tao)) - dt * exp(-dt / tao)) / t0;
	t3 = (1 - exp(-dt / tao)) / t0;
	t4 = (dt * exp(-dt / tao) - tao * (1 - exp(-dt / tao))) / t0;
	t5 = -tao * (1 - exp(-dt / tao)) / t0;

	//>spatial expansion of the equilibrium
	for (int i = 0; i < pnum; i++) {
		w_l[i] = (wl[i] - wcl[i]) / dil;
		w_r[i] = (wcr[i] - wr[i]) / dir;
		w_L[i] = (wo[i] - wcl[i]) / dil;
		w_R[i] = (wcr[i] - wo[i]) / dir;
	}
	get_slope(a_l, priml, w_l, 1);
	get_slope(a_r, primr, w_r, 1);
	get_slope(a_L, primo, w_L, 1);
	get_slope(a_R, primo, w_R, 1);

	//>time expansion of the equilibrium
	moment_a(mom_l, a_l, 1, 0, 0, mu_l, mv_l, mw_l, mi_l, priml);
	moment_a(mom_r, a_r, 1, 0, 0, mu_r, mv_r, mw_r, mi_r, primr);
	get_slope(A_l, priml, mom_l, -1);
	get_slope(A_r, primr, mom_r, -1);

	moment_uvw(mom_o, 0, 0, 0, 0, mu, mv, mw, mi, primo);
	moment_uvw(mom_l, 0, 0, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_uvw(mom_r, 0, 0, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);

	moment_a(mom_aL, a_L, 1, 0, 0, mul, mv, mw, mi, primo);
	moment_a(mom_aR, a_R, 1, 0, 0, mur, mv, mw, mi, primo);
	moment_a(mom_al, a_l, 1, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_a(mom_ar, a_r, 1, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);
	moment_a(mom_Al, A_l, 0, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_a(mom_Ar, A_r, 0, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);

	for (int i = 0; i < pnum; i++) {
		w_T[i] = t1 * mom_o[i] + t2 * (mom_aL[i] + mom_aR[i]) + t3 * (mom_l[i] + mom_r[i]) + (t4 + t5) * (mom_al[i] + mom_ar[i]) + t5 * (mom_Al[i] + mom_Ar[i]);
	}
	get_slope(A, primo, w_T, 1);

	//>flux for conservative variables
	t0 = dt + tao * exp(-dt / tao) - tao;
	t1 = -tao * dt - 2 * tao * tao * exp(-dt / tao) + 2 * tao * tao - tao * dt * exp(-dt / tao);
	t2 = 0.5 * dt * dt - dt * tao - tao * tao * exp(-dt / tao) + tao * tao;
	t3 = tao * (1 - exp(-dt / tao));
	t4 = tao * tao * exp(-dt / tao) - 2 * tao * tao + (tao * tao + tao * dt) * exp(-dt / tao);
	t5 = -tao * tao * (1 - exp(-dt / tao));

	moment_uvw(mom_o, 1, 0, 0, 0, mu, mv, mw, mi, primo);
	moment_uvw(mom_l, 1, 0, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_uvw(mom_r, 1, 0, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);

	moment_a(mom_aL, a_L, 2, 0, 0, mul, mv, mw, mi, primo);
	moment_a(mom_aR, a_R, 2, 0, 0, mur, mv, mw, mi, primo);
	moment_a(mom_al, a_l, 2, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_a(mom_ar, a_r, 2, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);
	moment_a(mom_Al, A_l, 1, 0, 0, mul_l, mv_l, mw_l, mi_l, priml);
	moment_a(mom_Ar, A_r, 1, 0, 0, mur_r, mv_r, mw_r, mi_r, primr);
	moment_a(mom_A, A, 1, 0, 0, mu, mv, mw, mi, primo);

	for (int i = 0; i < pnum; i++) {
		flux[i] = t0 * mom_o[i] + t1 * (mom_aL[i] + mom_aR[i]) + t2 * mom_A[i] + t3 * (mom_l[i] + mom_r[i]) + t4 * (mom_al[i] + mom_ar[i]) + t5 * (mom_Al[i] + mom_Ar[i]);
	}

	//>Pr fix
	if (Pr != 1) {
		double fde, ufde, vfde, wfde, q, q0, q1, q2, q3, u1, u2, u3;
		fde = primo[0] - tao * mompr(a_L, 1, 0, 0, 0, mu, mv, mw, mi, primo) + (dt - tao) * mompr(A, 0, 0, 0, 0, mu, mv, mw, mi, primo);
		ufde = primo[0] * primo[1] - tao * mompr(a_L, 2, 0, 0, 0, mu, mv, mw, mi, primo) + (dt - tao) * mompr(A, 1, 0, 0, 0, mu, mv, mw, mi, primo);
		vfde = primo[0] * primo[2] - tao * mompr(a_L, 1, 1, 0, 0, mu, mv, mw, mi, primo) + (dt - tao) * mompr(A, 0, 1, 0, 0, mu, mv, mw, mi, primo);
		wfde = primo[0] * primo[3] - tao * mompr(a_L, 1, 0, 1, 0, mu, mv, mw, mi, primo) + (dt - tao) * mompr(A, 0, 0, 1, 0, mu, mv, mw, mi, primo);

		q0 = mompr(a_L, 4, 0, 0, 0, mu, mv, mw, mi, primo) + mompr(A, 3, 0, 0, 0, mu, mv, mw, mi, primo) \
			+ mompr(a_L, 2, 2, 0, 0, mu, mv, mw, mi, primo) + mompr(A, 1, 2, 0, 0, mu, mv, mw, mi, primo) \
			+ mompr(a_L, 2, 0, 2, 0, mu, mv, mw, mi, primo) + mompr(A, 1, 0, 2, 0, mu, mv, mw, mi, primo) \
			+ mompr(a_L, 2, 0, 0, 2, mu, mv, mw, mi, primo) + mompr(A, 1, 0, 0, 2, mu, mv, mw, mi, primo);
		q1 = mompr(a_L, 3, 0, 0, 0, mu, mv, mw, mi, primo) + mompr(A, 2, 0, 0, 0, mu, mv, mw, mi, primo);
		q2 = mompr(a_L, 2, 1, 0, 0, mu, mv, mw, mi, primo) + mompr(A, 1, 1, 0, 0, mu, mv, mw, mi, primo);
		q3 = mompr(a_L, 2, 0, 1, 0, mu, mv, mw, mi, primo) + mompr(A, 1, 0, 1, 0, mu, mv, mw, mi, primo);

		u1 = ufde / fde; u2 = vfde / fde; u3 = wfde / fde;
		q = -tao * (q0 - u1 * q1 - u2 * q2 - u3 * q3) * dt;
		flux[4] = flux[4] + (1 / Pr - 1) * q;
	}

}

//--------------------------------------------------------------------------------
//>SOLVER
//--------------------------------------------------------------------------------
double dens_new[cnum], xmom_new[cnum], ymom_new[cnum], zmom_new[cnum], energy_new[cnum];
double dens_flux[snum * cnum], xmom_flux[snum * cnum], ymom_flux[snum * cnum], zmom_flux[snum * cnum], energy_flux[snum * cnum];

//>function: calculate time step
void get_time() {
	double v, sos;
	dt = 10.0;
	for (int i = 0; i < cnum; i++) {
		sos = (xmom[i * rnum] * xmom[i * rnum] + ymom[i * rnum] * ymom[i * rnum] + zmom[i * rnum] * zmom[i * rnum]) / dens[i * rnum];
		v = sqrt(sos / dens[i * rnum]);
		sos = 0.5 * dens[i * rnum] / (gam - 1.0) / (energy[i * rnum] - 0.5 * sos);
		sos = sqrt(0.5 * gam / sos);
		dt = min(sqrt(scale[(1 + snum) * i]) / (v + sos), dt);
	}
	dt = cfl * dt;
	if (sim_time + dt > max_time) { dt = max_time - sim_time; }
}

//>main function
void solver() {

	double flux[pnum], wcl[pnum], wcr[pnum], wl[pnum], wr[pnum], dil, dir, x[3], y[3], z[3];
	double area, volume, normal_1, normal_2, normal_3;
	double prim_cl[pnum], prim_cr[pnum];
	double res_dens = 0, res_xmom = 0, res_ymom = 0, res_zmom = 0, res_energy = 0;
	double avg_dens = 0, avg_xmom = 0, avg_ymom = 0, avg_zmom = 0, avg_energy = 0;
	int near, side_r;

	//>time step
	get_time();

	for (int i = 0; i < snum * cnum; i++) {
		dens_flux[i] = -10000;
		xmom_flux[i] = -10000;
		ymom_flux[i] = -10000;
		zmom_flux[i] = -10000;
		energy_flux[i] = -10000;
	}

	//>for each cell
	for (int i = 0; i < cnum; i++) {
		dens_new[i] = 0;
		xmom_new[i] = 0;
		ymom_new[i] = 0;
		zmom_new[i] = 0;
		energy_new[i] = 0;

		wcl[0] = dens[i * rnum];
		wcl[1] = xmom[i * rnum];
		wcl[2] = ymom[i * rnum];
		wcl[3] = zmom[i * rnum];
		wcl[4] = energy[i * rnum];

		volume = scale[(1 + snum) * i];
		x[0] = center[3 * i];
		y[0] = center[3 * i + 1];
		z[0] = center[3 * i + 2];

		//>for each side
		for (int j = 0; j < snum; j++) {

			wl[0] = dens[i * rnum + 1 + j];
			wl[1] = xmom[i * rnum + 1 + j];
			wl[2] = ymom[i * rnum + 1 + j];
			wl[3] = zmom[i * rnum + 1 + j];
			wl[4] = energy[i * rnum + 1 + j];

			area = scale[(1 + snum) * i + 1 + j];
			near = connect[i * snum + j];
			for (int k = 0; k < snum; k++) {
				if (near == -1) { side_r = -1; break; }
				if (connect[near * snum + k] == i) { side_r = k; break; }
			}

			x[1] = points[3 * snum * i + 3 * j];
			y[1] = points[3 * snum * i + 3 * j + 1];
			z[1] = points[3 * snum * i + 3 * j + 2];
			dil = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2) + pow(z[1] - z[0], 2));

			//>if boundary cells
			if (near == -1) {
				if (j == 0) {//----------------------------supersonic inlet
					prim_cr[0] = rho_in;
					prim_cr[1] = u_in;
					prim_cr[2] = v_in;
					prim_cr[3] = w_in;
					prim_cr[4] = 0.5 * rho_in / p_in;
					get_conserved(wcr, prim_cr);
					for (int k = 0; k < pnum; k++) { wr[k] = wcr[k]; }
				}
				else if (j == 1) {//----------------------supersonic outlet
					for (int k = 0; k < pnum; k++) { wcr[k] = wcl[k]; wr[k] = wcr[k]; }
				}
				else if (j == 2) {//-----------------------------lower wall
					if (center[3 * i] < 0.01) {
						wcr[0] = wcl[0];
						wcr[1] = wcl[1];
						wcr[2] = -wcl[2];
						wcr[3] = wcl[3];
						wcr[4] = wcl[4];
						for (int k = 0; k < pnum; k++) { wr[k] = wcr[k]; }
					}
					else {
						normal_1 = normal[3 * snum * i + 3 * j];
						normal_2 = normal[3 * snum * i + 3 * j + 1];
						normal_3 = normal[3 * snum * i + 3 * j + 2];
						rotation(wcl, normal_1, normal_2, normal_3, 1);
						wcr[0] = wcl[0];
						wcr[1] = -wcl[1];
						wcr[2] = wcl[2];
						wcr[3] = wcl[3];
						wcr[4] = wcl[4];
						rotation(wcl, normal_1, normal_2, normal_3, -1);
						rotation(wcr, normal_1, normal_2, normal_3, -1);
						for (int k = 0; k < pnum; k++) { wr[k] = wcr[k]; }
					}
				}
				else if (j == 3) {//-----------------------------upper wall
					if (center[3 * i] < 0.383868) {
						for (int k = 0; k < pnum; k++) { wcr[k] = wcl[k]; wr[k] = wcr[k]; }
					}
					else {
						wcr[0] = wcl[0];
						wcr[1] = wcl[1];
						wcr[2] = -wcl[2];
						wcr[3] = wcl[3];
						wcr[4] = wcl[4];
						for (int k = 0; k < pnum; k++) { wr[k] = wcr[k]; }
					}
				}
				else {//--------------------------------------side walls 
					wcr[0] = wcl[0];
					wcr[1] = wcl[1];
					wcr[2] = wcl[2];
					wcr[3] = -wcl[3];
					wcr[4] = wcl[4];
					for (int k = 0; k < pnum; k++) { wr[k] = wcr[k]; }
				}
				dir = dil;
				normal_1 = normal[3 * snum * i + 3 * j];
				normal_2 = normal[3 * snum * i + 3 * j + 1];
				normal_3 = normal[3 * snum * i + 3 * j + 2];

				rotation(wcl, normal_1, normal_2, normal_3, 1);
				rotation(wcr, normal_1, normal_2, normal_3, 1);
				rotation(wl, normal_1, normal_2, normal_3, 1);
				rotation(wr, normal_1, normal_2, normal_3, 1);
				fluxGKS(flux, wcl, wcr, wl, wr, dil, dir);
				rotation(flux, normal_1, normal_2, normal_3, -1);
				rotation(wcl, normal_1, normal_2, normal_3, -1);
			}
			//>if inner cells
			else {
				if (dens_flux[snum * i + j] == -10000) {
					wcr[0] = dens[near * rnum];
					wcr[1] = xmom[near * rnum];
					wcr[2] = ymom[near * rnum];
					wcr[3] = zmom[near * rnum];
					wcr[4] = energy[near * rnum];

					wr[0] = dens[near * rnum + 1 + side_r];
					wr[1] = xmom[near * rnum + 1 + side_r];
					wr[2] = ymom[near * rnum + 1 + side_r];
					wr[3] = zmom[near * rnum + 1 + side_r];
					wr[4] = energy[near * rnum + 1 + side_r];

					x[2] = center[3 * near];
					y[2] = center[3 * near + 1];
					z[2] = center[3 * near + 2];
					dir = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2) + pow(z[2] - z[1], 2));

					normal_1 = normal[3 * snum * i + 3 * j];
					normal_2 = normal[3 * snum * i + 3 * j + 1];
					normal_3 = normal[3 * snum * i + 3 * j + 2];

					rotation(wcl, normal_1, normal_2, normal_3, 1);
					rotation(wcr, normal_1, normal_2, normal_3, 1);
					rotation(wl, normal_1, normal_2, normal_3, 1);
					rotation(wr, normal_1, normal_2, normal_3, 1);
					fluxGKS(flux, wcl, wcr, wl, wr, dil, dir);
					rotation(flux, normal_1, normal_2, normal_3, -1);
					rotation(wcl, normal_1, normal_2, normal_3, -1);

					dens_flux[snum * i + j] = flux[0];
					xmom_flux[snum * i + j] = flux[1];
					ymom_flux[snum * i + j] = flux[2];
					zmom_flux[snum * i + j] = flux[3];
					energy_flux[snum * i + j] = flux[4];

					dens_flux[snum * near + side_r] = -1 * flux[0];
					xmom_flux[snum * near + side_r] = -1 * flux[1];
					ymom_flux[snum * near + side_r] = -1 * flux[2];
					zmom_flux[snum * near + side_r] = -1 * flux[3];
					energy_flux[snum * near + side_r] = -1 * flux[4];
				}
				else {
					flux[0] = dens_flux[snum * i + j];
					flux[1] = xmom_flux[snum * i + j];
					flux[2] = ymom_flux[snum * i + j];
					flux[3] = zmom_flux[snum * i + j];
					flux[4] = energy_flux[snum * i + j];
				}
			}

			//>volume conservation 
			dens_new[i] = dens_new[i] + flux[0] * area / volume;
			xmom_new[i] = xmom_new[i] + flux[1] * area / volume;
			ymom_new[i] = ymom_new[i] + flux[2] * area / volume;
			zmom_new[i] = zmom_new[i] + flux[3] * area / volume;
			energy_new[i] = energy_new[i] + flux[4] * area / volume;
		}
	}

	for (int i = 0; i < cnum; i++) {

		//>update cell-averaged variables
		dens[rnum * i] = dens[rnum * i] - dens_new[i];
		xmom[rnum * i] = xmom[rnum * i] - xmom_new[i];
		ymom[rnum * i] = ymom[rnum * i] - ymom_new[i];
		zmom[rnum * i] = zmom[rnum * i] - zmom_new[i];
		energy[rnum * i] = energy[rnum * i] - energy_new[i];

		//>record residual
		res_dens = res_dens + dens_new[i] * dens_new[i];
		res_xmom = res_xmom + xmom_new[i] * xmom_new[i];
		res_ymom = res_ymom + ymom_new[i] * ymom_new[i];
		res_zmom = res_zmom + zmom_new[i] * zmom_new[i];
		res_energy = res_energy + energy_new[i] * energy_new[i];

		avg_dens = avg_dens + fabs(dens[rnum * i]);
		avg_xmom = avg_xmom + fabs(xmom[rnum * i]);
		avg_ymom = avg_ymom + fabs(ymom[rnum * i]);
		avg_zmom = avg_zmom + fabs(zmom[rnum * i]);
		avg_energy = avg_energy + fabs(energy[rnum * i]);
	}
	res_dens = sqrt(cnum * res_dens) / avg_dens;
	res_xmom = sqrt(cnum * res_xmom) / avg_xmom;
	res_ymom = sqrt(cnum * res_ymom) / avg_ymom;
	res_zmom = sqrt(cnum * res_zmom) / avg_zmom;
	res_energy = sqrt(cnum * res_energy) / avg_energy;

	res = max(res_dens, res_energy);
	res = max(res, res_xmom);
	res = max(res, res_ymom);
	res = max(res, res_zmom);
}

//--------------------------------------------------------------------------------
//>INITIALIZATION
//-------------------------------------------------------------------------------- 
void initial() {

	//>set parameters
	cfl = 0.4;
	max_time = 5000;
	eps = 1.0E-5;
	output = 50;

	ck = 2;
	ct = 3;
	Pr = 1.0;
	gam = get_gamma();
	mu_ref = 1.187E-3;

	rho_in = 2.82;   u_in = 1725.6; v_in = 0.0;  w_in = 0.0;  p_in = 101325;
	rho_int = rho_in; u_int = u_in;   v_int = v_in; w_int = w_in; lam_int = 0.5 * rho_in / p_in;

	//>initialize the flow field
	double prim_int[pnum] = { rho_int,u_int,v_int,w_int,lam_int };
	double w_int[pnum];
	get_conserved(w_int, prim_int);

	for (int i = 0; i < cnum; i++) {
		dens[rnum * i] = w_int[0];
		xmom[rnum * i] = w_int[1];
		ymom[rnum * i] = w_int[2];
		zmom[rnum * i] = w_int[3];
		energy[rnum * i] = w_int[4];
	}

}

//--------------------------------------------------------------------------------
//>MAIN FUNCTION
//--------------------------------------------------------------------------------
int main() {

	read_mesh();
	no_zero();
	geometry_hex();

	initial();

	iter = 1;
	sim_time = 0;
	int write = 1;
	double sos, ma;

	while (sim_time < max_time) {

		reconstruct_1st();
		solver();
		sim_time = sim_time + dt;
		cout << "iter = " << iter << ", time = " << sim_time << ", res = " << res << endl;
		iter++;

		//-------------------------------------------
		if (write == output) {
			ofstream outfile("mesh/hypersonic_intake/out_ma.dat");
			outfile << "variables=\"x\",\"y\",\"z\",\"ma\"" << endl;
			outfile << "zone N=" << nnum << ",E=" << cnum << ", datapacking=block" << endl;
			outfile << "varlocation=([4]=cellcentered)" << endl;
			outfile << "zonetype=FEBRICK" << endl;
			for (int i = 0; i < nnum; i++) { outfile << notes[3 * i] << endl; }
			for (int i = 0; i < nnum; i++) { outfile << notes[3 * i + 1] << endl; }
			for (int i = 0; i < nnum; i++) { outfile << notes[3 * i + 2] << endl; }
			for (int i = 0; i < cnum; i++) {
				ma = (pow(xmom[i * rnum], 2) + pow(ymom[i * rnum], 2) + pow(zmom[i * rnum], 2)) / dens[i * rnum];
				sos = 0.5 * dens[i * rnum] / (gam - 1.0) / (energy[i * rnum] - 0.5 * ma);
				sos = sqrt(0.5 * gam / sos);
				outfile << sqrt(ma / dens[i * rnum]) / sos << endl;
			}
			for (int i = 0; i < cnum; i++) { outfile << element[i * vnum] + 1 << " " << element[i * vnum + 1] + 1 << " " << element[i * vnum + 2] + 1 << " " << element[i * vnum + 3] + 1 << " " << element[i * vnum + 4] + 1 << " " << element[i * vnum + 5] + 1 << " " << element[i * vnum + 6] + 1 << " " << element[i * vnum + 7] + 1 << endl; }
			outfile.close();

			ofstream outfile1("mesh/hypersonic_intake/out_cp.dat");
			outfile1 << "variables=\"x\",\"y\",\"z\",\"cp\"" << endl;
			outfile1 << "zone N=" << nnum << ",E=" << cnum << ", datapacking=block" << endl;
			outfile1 << "varlocation=([4]=cellcentered)" << endl;
			outfile1 << "zonetype=FEBRICK" << endl;
			for (int i = 0; i < nnum; i++) { outfile1 << notes[3 * i] << endl; }
			for (int i = 0; i < nnum; i++) { outfile1 << notes[3 * i + 1] << endl; }
			for (int i = 0; i < nnum; i++) { outfile1 << notes[3 * i + 2] << endl; }
			for (int i = 0; i < cnum; i++) { outfile1 << ((gam - 1.0) * (energy[i * rnum] - 0.5 * (xmom[i * rnum] * xmom[i * rnum] + ymom[i * rnum] * ymom[i * rnum] + zmom[i * rnum] * zmom[i * rnum]) / dens[i * rnum]) - p_in) * 2 / rho_in / u_in / u_in << endl; }
			for (int i = 0; i < cnum; i++) { outfile1 << element[i * vnum] + 1 << " " << element[i * vnum + 1] + 1 << " " << element[i * vnum + 2] + 1 << " " << element[i * vnum + 3] + 1 << " " << element[i * vnum + 4] + 1 << " " << element[i * vnum + 5] + 1 << " " << element[i * vnum + 6] + 1 << " " << element[i * vnum + 7] + 1 << endl; }
			outfile1.close();

			std::cout << "\nResults have been written.\n\n";

			write = 0;
		}
		write++;
	}

	return 0;
}
