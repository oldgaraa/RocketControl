// Интегрирование математической модели движения ракеты, расчет и реализация наведения для попадания в заданную цель
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include "math.h"
#include "GOST4401-81.h"
using namespace std;
ofstream txtfile2("2.txt");

double a[] = { 0.0,0.5,0.5,1.0 };
#define dtdt  K[0][j]
#define dxgdt  K[1][j]
#define dygdt  K[2][j]
#define dzgdt  K[3][j]
#define dVxgdt K[4][j]
#define dVygdt K[5][j]
#define dVzgdt K[6][j]
#define domegaxdt K[7][j]
#define domegaydt K[8][j]
#define domegazdt K[9][j]
#define droRGdt K[10][j]
#define dlamdaRGdt K[11][j]
#define dmuRGdt K[12][j]
#define dnuRGdt K[13][j]

#define t  par[0]
#define xg  par[1]
#define yg  par[2]
#define zg  par[3]
#define Vxg par[4]
#define Vyg par[5]
#define Vzg par[6]
#define omegax par[7]
#define omegay par[8]
#define omegaz par[9]
#define roRG par[10]
#define lamdaRG par[11]
#define muRG par[12]
#define nuRG par[13]

#define t_RK  parbuf[0]
#define xg_RK  parbuf[1]
#define yg_RK  parbuf[2]
#define zg_RK  parbuf[3]
#define Vxg_RK parbuf[4]
#define Vyg_RK parbuf[5]
#define Vzg_RK parbuf[6]
#define omegax_RK parbuf[7]
#define omegay_RK parbuf[8]
#define omegaz_RK parbuf[9]
#define roRG_RK parbuf[10]
#define lamdaRG_RK parbuf[11]
#define muRG_RK parbuf[12]
#define nuRG_RK parbuf[13]


// Объявление функций
double Cx(double M);
double Cy_alpha(double mach);
double Cz_betta(double M, double betta);
double mz_alpha(double M);
double my_betta(double M);
void MATR(double ro, double lamda, double mu, double nu, double(*A)[3][3], double(*AT)[3][3]);
void Integrator(double(*parInteg)[14], const double h);
const double g = 9.81;
double K[99][99], par[99], parbuf[99], par_drob[99],
A[3][3] = { 0 }, AT[3][3] = { 0 }, B[3][3] = { 0 };
double tang, risk, kren, m, V0, q, M, r, V, X, Y, Z, Mx, My, Mz, Mstab, alpha, betta, alpha_full, vx, vy, vz, vr, vfi, vhi,
delta_tang, delta_risk, delta_kren, fi, hi, dfidt, dhidt, x_t, y_t, z_t, Vx_t = 0, Vy_t = 0, Vz_t = 0,
Ix = 14.8E-6,
Iy = 27E-6,
Iz = 27E-6,
dm = 9.53E-3,
Sm = M_PI*dm*dm / 4.0,
L = 0.00767,
dt;

double par0[14] = { 0 }; //вектор начальных параметров
int main()
{


	Atm_GOST();

			//Начальные условия
	tang = 1 * M_PI / 180.0;
			risk = 0;
			kren = 0;
			m = 3.2E-3;
			V0 = 228;
			dt = 0.0001;
			t = 0;
			xg = 0;
			yg = 1.5;
			zg = 0;
			Vxg = V0*cos(tang);
			Vyg = V0*sin(tang);
			Vzg = 0;
			omegax = 2 * M_PI*V0 / 0.19; //2 * M_PI*V0 / (52 / 6);
			omegay = 0;
			omegaz = 0;
			roRG = cos(risk*0.5)*cos(tang*0.5)*cos(kren*0.5) - sin(risk*0.5)*sin(kren*0.5)*sin(tang*0.5);
			lamdaRG = sin(risk*0.5)*sin(tang*0.5)*cos(kren*0.5) + cos(risk*0.5)*cos(tang*0.5)*sin(kren*0.5);
			muRG = sin(risk*0.5)*cos(tang*0.5)*cos(kren*0.5) + cos(risk*0.5)*sin(tang*0.5)*sin(kren*0.5);
			nuRG = cos(risk*0.5)*sin(tang*0.5)*cos(kren*0.5) - sin(risk*0.5)*cos(tang*0.5)*sin(kren*0.5);
			//for (int i : par) par0[i] = par[i]; //фиксация начальных параметров

			ofstream txtfile("1.txt");

			txtfile << "t,s" << ";" << "x,m" << ";" << "y,m" << ";" << "z,m" << ";" << "Vx,m/s" << ";" << "Vy,m/s" << ";" << "Vz,m/s" << ";" << "V,m/s" << ";"
				<< "tang,grad" << ";" << "risk,grad" << ";" << "kren,grad" << ";" << "omega_x,grad/s" << ";" << "omega_y,grad/s" << ";" << "omega_z,grad/s" << ";"
				<< "alpha,grad" << ";" << "betta,grad" << ";" << "deltaT, drad" << ";" << "deltaN, grag" << ";" << "deltaE, grad" << ";" << endl;

			// Интегрирование
			int step = 0;
			for (;;) {
				// Расчет промежуточных параметров на шаге
				MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
				tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
				risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
				kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
				V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
					q = fro(yg)*V*V*0.5,
					vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
					vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
					vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
					alpha = -atan2(vy, vx);
				betta = asin(vz / V);
				alpha_full = pow(alpha*alpha + betta*betta, 0.5);
				M = V / fa(yg);
				X = -Cx(M)*q*Sm;
				Y = (Cy_alpha(M)*alpha)*q*Sm;
				Z = (Cz_betta(M, betta)*betta)*q*Sm;
				Mx = 0;
				My = (my_betta(M)*betta)*q*Sm*L;
				Mz = (mz_alpha(M)*alpha)*q*Sm*L;

				// Дробление шага
				if (yg < 0) {
					dt *= 0.5;
					for (int i = 0; i <= 13; i++) par[i] = par_drob[i];
				}
				else {
					for (int i = 0; i <= 13; i++) par_drob[i] = par[i];
					//вывод в файл
					int st = 0.0001 / dt;
					if (step % st == 0) 
					txtfile
						<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
						<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
						<< omegax << ";" << omegay << ";" << omegaz << ";"
						<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
						<< delta_tang * 180 / M_PI << ";" << delta_risk * 180 / M_PI << ";" << delta_kren * 180 / M_PI << ";" << r << ";" << endl;
					cout << yg << endl;
				}
				if (abs(yg) < 0.0001) break;	//Условие выхода
				step++;
				//Цикл расчета коэффициентов
				for (int j = 1; j <= 4; j++) {
					for (int i = 0; i <= 13; i++)  parbuf[i] = par[i] + K[i][j - 1] * a[j - 1]; // Считаем новые аргументы

																								// Нормировка параметров РГ
					double norm_RG = pow(roRG_RK*roRG_RK + lamdaRG_RK*lamdaRG_RK + nuRG_RK*nuRG_RK + muRG_RK*muRG_RK, 0.5);
					roRG_RK /= norm_RG;
					lamdaRG_RK /= norm_RG;
					nuRG_RK /= norm_RG;
					muRG_RK /= norm_RG;
					MATR(roRG_RK, lamdaRG_RK, muRG_RK, nuRG_RK, &A, &AT);

					// Расчет промежуточных параметров на шаге
					MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
					tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
					risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
					kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
					V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
						q = fro(yg)*V*V*0.5,
						vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
						vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
						vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
						alpha = -atan2(vy, vx);
					betta = asin(vz / V);
					alpha_full = pow(alpha*alpha + betta*betta, 0.5);
					M = V / fa(yg);
					X = -Cx(M)*q*Sm;
					Y = (Cy_alpha(M)*alpha)*q*Sm;
					Z = (Cz_betta(M, betta)*betta)*q*Sm;
					Mx = 0;
					My = (my_betta(M)*betta)*q*Sm*L;
					Mz = (mz_alpha(M)*alpha)*q*Sm*L;

					//Производные (коэффициенты K[i][j])
					dtdt = 1 * dt,
						dxgdt = Vxg_RK*dt,
						dygdt = Vyg_RK*dt,
						dzgdt = Vzg_RK*dt,
						dVxgdt = 1 / m*(X*A[0][0] + Y*A[0][1] + Z*A[0][2])*dt,
						dVygdt = 1 / m*(X*A[1][0] + Y*A[1][1] + Z*A[1][2] - m*g)*dt,
						dVzgdt = 1 / m*(X*A[2][0] + Y*A[2][1] + Z*A[2][2])*dt,
						domegaxdt = (Mx / Ix - omegay_RK*omegaz_RK*(Iz - Iy) / Ix)*dt,
						domegaydt = (My / Iy - omegax_RK*omegaz_RK*(Ix - Iz) / Iy)*dt,
						domegazdt = (Mz / Iz - omegay_RK*omegax_RK*(Iy - Ix) / Iz)*dt,
						droRGdt = -(omegax_RK*lamdaRG_RK + omegay_RK*muRG_RK + omegaz_RK*nuRG_RK)*0.5*dt,
						dlamdaRGdt = (omegax_RK*roRG_RK - omegay_RK*nuRG_RK + omegaz_RK*muRG_RK)*0.5*dt,
						dmuRGdt = (omegax_RK*nuRG_RK + omegay_RK*roRG_RK - omegaz_RK*lamdaRG_RK)*0.5*dt,
						dnuRGdt = (-omegax_RK*muRG_RK + omegay_RK*lamdaRG_RK + omegaz_RK*roRG_RK)*0.5*dt;

				}

				// Приращения
				for (int i = 0; i <= 13; i++) {
					double b[5] = { 0.0, 1.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0 };
					for (int j = 1; j <= 4; j++) par[i] += K[i][j] * b[j];
				}
				// Нормировка параметров РГ
				double norm_RG = pow(roRG*roRG + lamdaRG*lamdaRG + nuRG*nuRG + muRG*muRG, 0.5);
				roRG /= norm_RG;
				lamdaRG /= norm_RG;
				nuRG /= norm_RG;
				muRG /= norm_RG;
			}  // Конец интегрирования траектории
			   // Вывод последнего шага
			txtfile
				<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
				<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
				<< omegax << ";" << omegay << ";" << omegaz << ";"
				<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
				<< delta_tang * 180 / M_PI << ";" << delta_risk * 180 / M_PI << ";" << delta_kren * 180 / M_PI << ";" << r << ";" << endl;


			//Integrator(par, dt);



	cout << "END!" << endl;
	system("pause");
    return 0;
}

//Аэродинамика
double Cx(double M) {
	return 0.20;
}
double Cy_alpha(double M) {
	return 0.70 / (5 * M_PI / 180.0);
}
double Cz_betta(double M, double betta) {
	return -Cy_alpha(M);
}
double mz_alpha(double M) {
	return 3E-4/ (5 * M_PI / 180.0);
}
double my_betta(double M) {
	return mz_alpha(M);
}

//Матрицы
void MATR(double ro, double lamda, double mu, double nu, double(*A)[3][3], double(*AT)[3][3]) {
	double a[3][3] = { 0 };
	// Матрица из нормальной земной в связную
	a[0][0] = ro*ro + lamda*lamda - mu*mu - nu*nu;
	a[0][1] = 2 * (-ro*nu + lamda*mu);
	a[0][2] = 2 * (ro*mu + lamda*nu);

	a[1][0] = 2 * (ro*nu + lamda*mu);
	a[1][1] = ro*ro - lamda*lamda + mu*mu - nu*nu;
	a[1][2] = 2 * (-ro*lamda + nu*mu);

	a[2][0] = 2 * (-ro*mu + lamda*nu);
	a[2][1] = 2 * (ro*lamda + nu*mu);
	a[2][2] = ro*ro - lamda*lamda - mu*mu + nu*nu;

	for (int i = 0; i <= 2; i++) for (int j = 0; j <= 2; j++)
		(*A)[i][j] = a[i][j];
	//Из связной в нормальную
	for (int i = 0; i <= 2; i++) for (int j = 0; j <= 2; j++)
		(*AT)[i][j] = a[j][i];
}

void Integrator(double (*parInteg)[13], const double h) {
	for (int i = 0; i <= 13; i++) par[i] = (*parInteg)[i];
	double dt = h;

	ofstream txtfile("1.txt");

	txtfile << "t,s" << ";" << "x,m" << ";" << "y,m" << ";" << "z,m" << ";" << "Vx,m/s" << ";" << "Vy,m/s" << ";" << "Vz,m/s" << ";" << "V,m/s" << ";"
		<< "tang,grad" << ";" << "risk,grad" << ";" << "kren,grad" << ";" << "omega_x,grad/s" << ";" << "omega_y,grad/s" << ";" << "omega_z,grad/s" << ";"
		<< "alpha,grad" << ";" << "betta,grad" << ";" << "deltaT, drad" << ";" << "deltaN, grag" << ";" << "deltaE, grad" << ";" << endl;

	// Интегрирование
	int step = 0;
	for (;;) {
		// Расчет промежуточных параметров на шаге
		MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
		tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
		risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
		kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
		V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
			q = fro(yg)*V*V*0.5,
			vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
			vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
			vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
			alpha = -atan2(vy, vx);
		betta = asin(vz / V);
		alpha_full = pow(alpha*alpha + betta*betta, 0.5);
		M = V / fa(yg);
		X = -Cx(M)*q*Sm;
		Y = (Cy_alpha(M)*alpha)*q*Sm;
		Z = (Cz_betta(M, betta)*betta)*q*Sm;
		Mx = 0;
		My = (my_betta(M)*betta)*q*Sm*L;
		Mz = (mz_alpha(M)*alpha)*q*Sm*L;

		// Дробление шага
		if (yg < 0) {
			dt *= 0.5;
			for (int i = 0; i <= 13; i++) par[i] = par_drob[i];
		}
		else {
			for (int i = 0; i <= 13; i++) par_drob[i] = par[i];
			//вывод в файл
			int st = 1 / dt;
			//if (step % st == 0) 
			txtfile
				<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
				<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
				<< omegax << ";" << omegay << ";" << omegaz << ";"
				<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
				<< delta_tang * 180 / M_PI << ";" << delta_risk * 180 / M_PI << ";" << delta_kren * 180 / M_PI << ";" << r << ";" << endl;
		}
		if (abs(yg) < 0.0001) break;	//Условие выхода
		step++;
		//Цикл расчета коэффициентов
		for (int j = 1; j <= 4; j++) {
			for (int i = 0; i <= 13; i++)  parbuf[i] = par[i] + K[i][j - 1] * a[j - 1]; // Считаем новые аргументы

																						// Нормировка параметров РГ
			double norm_RG = pow(roRG_RK*roRG_RK + lamdaRG_RK*lamdaRG_RK + nuRG_RK*nuRG_RK + muRG_RK*muRG_RK, 0.5);
			roRG_RK /= norm_RG;
			lamdaRG_RK /= norm_RG;
			nuRG_RK /= norm_RG;
			muRG_RK /= norm_RG;
			MATR(roRG_RK, lamdaRG_RK, muRG_RK, nuRG_RK, &A, &AT);

			// Расчет промежуточных параметров на шаге
			MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
			tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
			risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
			kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
			V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
				q = fro(yg)*V*V*0.5,
				vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
				vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
				vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
				alpha = -atan2(vy, vx);
			betta = asin(vz / V);
			alpha_full = pow(alpha*alpha + betta*betta, 0.5);
			M = V / fa(yg);
			X = -Cx(M)*q*Sm;
			Y = (Cy_alpha(M)*alpha)*q*Sm;
			Z = (Cz_betta(M, betta)*betta)*q*Sm;
			Mx = 0;
			My = (my_betta(M)*betta)*q*Sm*L;
			Mz = (mz_alpha(M)*alpha)*q*Sm*L;

			//Производные (коэффициенты K[i][j])
			dtdt = 1 * dt,
				dxgdt = Vxg_RK*dt,
				dygdt = Vyg_RK*dt,
				dzgdt = Vzg_RK*dt,
				dVxgdt = 1 / m*(X*A[0][0] + Y*A[0][1] + Z*A[0][2])*dt,
				dVygdt = 1 / m*(X*A[1][0] + Y*A[1][1] + Z*A[1][2] - m*g)*dt,
				dVzgdt = 1 / m*(X*A[2][0] + Y*A[2][1] + Z*A[2][2])*dt,
				domegaxdt = (Mx / Ix - omegay_RK*omegaz_RK*(Iz - Iy) / Ix)*dt,
				domegaydt = (My / Iy - omegax_RK*omegaz_RK*(Ix - Iz) / Iy)*dt,
				domegazdt = (Mz / Iz - omegay_RK*omegax_RK*(Iy - Ix) / Iz)*dt,
				droRGdt = -(omegax_RK*lamdaRG_RK + omegay_RK*muRG_RK + omegaz_RK*nuRG_RK)*0.5*dt,
				dlamdaRGdt = (omegax_RK*roRG_RK - omegay_RK*nuRG_RK + omegaz_RK*muRG_RK)*0.5*dt,
				dmuRGdt = (omegax_RK*nuRG_RK + omegay_RK*roRG_RK - omegaz_RK*lamdaRG_RK)*0.5*dt,
				dnuRGdt = (-omegax_RK*muRG_RK + omegay_RK*lamdaRG_RK + omegaz_RK*roRG_RK)*0.5*dt;

		}

		// Приращения
		for (int i = 0; i <= 13; i++) {
			double b[5] = { 0.0, 1.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0 };
			for (int j = 1; j <= 4; j++) par[i] += K[i][j] * b[j];
		}
		// Нормировка параметров РГ
		double norm_RG = pow(roRG*roRG + lamdaRG*lamdaRG + nuRG*nuRG + muRG*muRG, 0.5);
		roRG /= norm_RG;
		lamdaRG /= norm_RG;
		nuRG /= norm_RG;
		muRG /= norm_RG;
	}  // Конец интегрирования траектории
	   // Вывод последнего шага
	txtfile
		<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
		<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
		<< omegax << ";" << omegay << ";" << omegaz << ";"
		<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
		<< delta_tang * 180 / M_PI << ";" << delta_risk * 180 / M_PI << ";" << delta_kren * 180 / M_PI << ";" << r << ";" << endl;

	for (int i = 0; i <= 13; i++) (*parInteg)[i] = par[i];
}