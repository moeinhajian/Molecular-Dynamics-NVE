#include <iostream>
#include <fstream>
#include <cmath>
#include<cstdlib>
#include<ctime>

using namespace std;
const int L = 12;                                           //length of cell
const int N = 1348;                                         // Number of atoms
const int mdts = 50000;                                     // Number of time steps
const double dt = 0.001;                                  // Time step
const double dr = 0.1;                                    // delta r for RDF
double T = 1.25;                                           // Temparature                                              // initial time
int x, y, z, i, j, b, a, w, q, d, cnt, nop;
double vi, cf, rtem, sumr=0, sumrsq;
double p[N][3], v[N][3], at[N][3], atdt[N][3];
double rx[N][N] = { 0 }, ry[N][N] = { 0 }, rz[N][N] = { 0 }, r[N][N] = { 0 }, g[60] = { 0 }, ginitial[60];
int main()
{
	srand((unsigned)time(NULL));
	ofstream outfile;

	// Generate regular initial location of each particle in the Cell
	// 12 atoms in X direction and 11 atoms in each Y ans Z directions until all 1348 matoms are displaced in the Cell with length 12 and the density equals to 0.78
	nop = 0;
	for (x = 0; x < 11; ++x) {
		for (y = 0; y < 11; ++y) {
			for (z = 0; z < 12; ++z) {
				if (cnt < N) {
					p[nop][0] = x * 1.091;
					p[nop][1] = y * 1.091;
					p[nop][2] = z * 1;
					++nop;
				}
			}
		}
	}
    double ktd = 0, vx = 0, vy = 0, vz = 0;
	// Velocity generation for every atoms
	for (i = 0; i < N; ++i) {
		v[i][0] = ((double)rand() / RAND_MAX) - 0.5;
		v[i][1] = ((double)rand() / RAND_MAX) - 0.5;
		v[i][2] = ((double)rand() / RAND_MAX) - 0.5;
		vi = sqrt(pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2));
		ktd += pow(vi, 2);
        // cout << "ktd=    " << vi << endl;
	}
    // cout << "ktd=    " << ktd << endl;
    // correction of initial velocity
    cf = sqrt(((double)3 * N * T) / ktd);
    ktd = 0;
    for (i = 0; i < N; ++i) {
        v[i][0] = cf * (v[i][0]);
        v[i][1] = cf * (v[i][1]);
        v[i][2] = cf * (v[i][2]);
        vi = sqrt(pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2));
        ktd += pow(vi, 2);
    }
    double uE = 0, kE, Et, E0;
	// computing distance between all molecules and potential energy of the Cell 
	for (i = 0; i < N - 1; ++i) {
		for (j = i + 1; j < N; ++j) {
			rx[i][j] = p[j][0] - p[i][0];
			ry[i][j] = p[j][1] - p[i][1];
			rz[i][j] = p[j][2] - p[i][2];
			// periodic boundary condition
			rx[i][j] -= (double)L * lround(rx[i][j] / L);
			ry[i][j] -= (double)L * lround(ry[i][j] / L);
			rz[i][j] -= (double)L * lround(rz[i][j] / L);

			r[i][j] = sqrt(pow(rx[i][j], 2) + pow(ry[i][j], 2) + pow(rz[i][j], 2));
			rx[j][i] = -rx[i][j];
			ry[j][i] = -ry[i][j];
			rz[j][i] = -rz[i][j];
			r[j][i] = r[i][j];
			if (r[i][j] <= 2.5) {
				uE += 4 * ((1 / pow((r[i][j]), 12) - 1 / pow((r[i][j]), 6)) - (1 / pow(2.5, 12) - 1 / pow(2.5, 6)));
			}
		}
	}
	
    cout << "ue=    " << uE << endl;
    
    
    uE = 0;
	// Calculating initial acceleration and put it at old acceleration array
	for (i = 0; i < N; ++i) {
		double Atem[3] = { 0 };
		for (j = 0; j < N; ++j) {
			if (j != i && r[i][j] <= 2.5) {
				Atem[0] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-rx[i][j] / r[i][j]);
				Atem[1] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-ry[i][j] / r[i][j]);
				Atem[2] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-rz[i][j] / r[i][j]);
                uE += 2 * ((1 / pow((r[i][j]), 12) - 1 / pow((r[i][j]), 6)) - (1 / pow(2.5, 12) - 1 / pow(2.5, 6)));
			}
		}
		at[i][0] = Atem[0];
		at[i][1] = Atem[1];
		at[i][2] = Atem[2];
	}

    
	outfile.open("Final1 energy-RDF.txt");
    kE = ktd / 2;
	E0 = uE + kE;
	d = 0; cnt = 0;
    double t = 0; 
	//main loop for all time steps
	while (t <= (dt * mdts)) {
		kE = 0, uE = 0, ktd = 0;
		// updating pinates according to the velocity verlet algorithm
		for (i = 0; i < N; ++i) {
			p[i][0] += dt * v[i][0] + pow(dt, 2) * at[i][0] / 2;
			p[i][1] += dt * v[i][1] + pow(dt, 2) * at[i][1] / 2;
			p[i][2] += dt * v[i][2] + pow(dt, 2) * at[i][2] / 2;
		}
		//calculating new distance between atoms and the potential energy for the system        
        for (i = 0; i < N - 1; ++i) {
			for (j = i + 1; j < N; ++j) {
				rx[i][j] = p[j][0] - p[i][0];
				ry[i][j] = p[j][1] - p[i][1];
				rz[i][j] = p[j][2] - p[i][2];
				rx[i][j] -= (double)L * lround(rx[i][j] / L);
				ry[i][j] -= (double)L * lround(ry[i][j] / L);
				rz[i][j] -= (double)L * lround(rz[i][j] / L);

				r[i][j] = sqrt(pow(rx[i][j], 2) + pow(ry[i][j], 2) + pow(rz[i][j], 2));
				rx[j][i] = -rx[i][j];
				ry[j][i] = -ry[i][j];
				rz[j][i] = -rz[i][j];
				r[j][i] = r[i][j];
				if (r[i][j] <= 2.5) {
					uE += 4 * ((1 / pow((r[i][j]), 12) - 1 / pow((r[i][j]), 6)) - (1 / pow(2.5, 12) - 1 / pow(2.5, 6)));
				}
			}
		}
        // cout << uE << endl;
		// calculating RDF in several time steps
		if ((d>500) && (d % 10 == 0)) {
			++cnt;
			for (i = 1; i <= 60; ++i) {
                q = 0;
                for (a = 0; a < N - 1; ++a) {
                    for (b = a + 1; b < N; ++b) {
                        rtem = r[a][b];
						sumr += rtem * rtem;
                        if (rtem <= (i * dr) && rtem > (((__int64)i - 1) * dr)) {
                            q = q + 1;
                        }
                    }
                }
                g[i - 1] += (__int64)q * 3 * pow(L, 3) / (2 * 3.14 * N * ((__int64)N - 1) * (pow((i * dr), 3) - pow((((__int64)i - 1) * dr), 3)));
            }
		}
		x = 0;
		if ((d % 10 == 0)) {
			for (a = 0; a < N - 1; ++a) {
				for (b = a + 1; b < N; ++b) {
					++x;
					rtem = r[a][b];
					sumr += rtem * rtem;
				}
			}
			sumrsq = sumr / x;
			cout<<"time: "<< t << " r mean: "<< sumrsq <<endl;
		}
		// Computing the new acceleration
		for (i = 0; i < N; ++i) {
			double Atem[3] = { 0 };
			for (j = 0; j < N; ++j) {
				if (j != i && r[i][j] <= 2.5) {
					Atem[0] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-rx[i][j] / r[i][j]);
					Atem[1] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-ry[i][j] / r[i][j]);
					Atem[2] += (24 / r[i][j]) * (2 * pow((1 / r[i][j]), 12) - pow((1 / r[i][j]), 6)) * (-rz[i][j] / r[i][j]);
				}
			}
			atdt[i][0] = Atem[0];
			atdt[i][1] = Atem[1];
			atdt[i][2] = Atem[2];
		}
        // cout << atdt << endl;
		// updating the velocity by velocity verlet for all atoms
		vx = 0;
		vy = 0;
		vz = 0;
		for (i = 0; i < N; ++i) {
			v[i][0] = v[i][0] + dt * (at[i][0] + atdt[i][0]) / 2;
			v[i][1] = v[i][1] + dt * (at[i][1] + atdt[i][1]) / 2;
			v[i][2] = v[i][2] + dt * (at[i][2] + atdt[i][2]) / 2;
            vi = sqrt(pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2));
			ktd = ktd + pow(vi, 2);
		}
		//rescaling velocity in fist 1000 steps
		if (d < 1000) {
			cf = sqrt(((double)3 * N * T) / ktd);
			ktd = 0;
			for (i = 0; i < N; ++i) {
				v[i][0] = v[i][0] * cf;
				v[i][1] = v[i][1] * cf;
				v[i][2] = v[i][2] * cf;
				vi = sqrt(pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2));
				ktd = ktd + pow(vi, 2);
			}
		}
		// computing kinetic energy and Temperature
		kE = ktd / 2;
		T = ktd / (3 * (__int64)N);
		//Replacing old acceleration by new acceleration
        // memcpy(at, atdt, sizeof(at));
		for (i = 0; i < N; ++i) {
			at[i][0] = atdt[i][0];
			at[i][1] = atdt[i][1];
			at[i][2] = atdt[i][2];
		}
		// computing the total energy
		Et = uE + kE;
		// printing the energies and temperatures in a text file
		outfile << " at time  " << t << "    Kt:   " << kE << "      Ut:   " << uE << "     E   " << Et << "     T-reduced  " << T << endl;
		t = t + dt;
		d += 1;
	}
	// printing final average RDF in the text file
	for (a = 0; a < 60; ++a) {
		g[a] = g[a] / cnt;
		outfile << "  between radius    " << ((a)*dr) << "  and    " << (((__int64)a + 1) * dr) << "  #particle    " << g[a] << endl;
	}

	outfile.close();

	return 0;
}
