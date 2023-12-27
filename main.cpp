#include <cmath>
#include <iostream>

using namespace std;

const int SIZE = 4;
const double A = -100;
const double B = -100;

// Function  k1 * x^3 + k2 * x^2 + k3 * x + k4
double f(double k[SIZE], double x)
{
	return k[0] * pow(x, 3) + k[1] * pow(x, 2) + k[2] * x + k[3];
}

// function integral k1 * x^3 + k2 * x^2 + k3 * x + k4
// Result (excluding constant)
//  k1/4 * x^4 + k2/3 * x^3 + k3/2 * x^2 + k4 * x

double integral(double k[SIZE], double x)
{
	return k[0] * pow(x, 4) * 0.25 + k[1] * pow(x, 3) / 3 + k[2] * pow(x, 2) * 0.5 + k[3] * x;
}

// defined integral from a to b
// Computed using the Newton-Leibniz formula
double definite_integral(double k[SIZE], double a, double b)
{
	return integral(k, b) - integral(k, a);
}

// distance between 2 function's grapphs on [a, b]
double dist_graph(double g1[SIZE], double g2[SIZE], double a, double b)
{
	double sum = 0;
	for (double x = a; x <= b; x += 0.01)
		sum += abs(f(g1, x) - f(g2, x));
	return sum;
}

// Distance between generations
double dist_gen(double g1[SIZE], double g2[SIZE])
{
	double d = 0;
	for (int i = 0; i < SIZE; i++)
		d += abs(g1[i] - g2[i]);
	return d;
}

// random existing  number from dMin to dMax
double fRand(double dMin, double dMax)
{
	double d = (double)rand() / RAND_MAX;
	return dMin + d * (dMax - dMin);
}

// New generation generating
void generate(double etalon[SIZE], double k[SIZE], double d, int N, int type, double new_k[4])
{
	for (int i = 0; i < SIZE; i++)
		new_k[i] = k[i] + fRand(-d, d);

	double dist;

	// Indirect sign
	if (!type)
		dist = abs(definite_integral(etalon, A, B)
			- definite_integral(new_k, A, B));
	// Main feature
	else
		dist = dist_graph(etalon, k, A, B);

	for (int i = 1; i < N; i++)
	{
		double new_k2[SIZE];
		for (int i = 0; i < SIZE; i++)
			new_k2[i] = k[i] + fRand(-d, d);
		double dist2;
		if (!type)
			dist2 = abs(definite_integral(new_k, A, B)
				- definite_integral(new_k2, A, B));
		else
			dist2 = dist_graph(new_k, new_k2, A, B);
		if (dist2 < dist)
		{
			dist = dist2;
			for (int i = 0; i < SIZE; i++)
				new_k[i] = new_k2[i];
		}
	}
}

// Running the evolutionary algorithm
double* runEvolutional(double et_k[SIZE], double k[SIZE], double d, int N, int type)
{
	double g1[SIZE], g2[SIZE];
	generate(et_k, k, d, N, type, g1);
	int num = 1;
	cout << "Generation " << num << ": ";
	for (int i = 0; i < SIZE; i++)
		cout << g1[i] << " ";
	cout << endl;
	while (true)
	{
		generate(et_k, g1, d, N, type, g2);
		cout << "Generation " << ++num << ": ";
		for (int i = 0; i < SIZE; i++)
			cout << g2[i] << " ";
		cout << endl;
		if (dist_gen(g1, g2) < 0.4 * d)
			break;
		for (int i = 0; i < SIZE; i++)
			g1[i] = g2[i];
	}
	return g2;
}

int main()
{

	// reference function
	double etalon[SIZE];
	etalon[0] = 1;
	etalon[1] = 2;
	etalon[2] = 3;
	etalon[3] = 4;

	double k[SIZE], d;
	int N;
	int type;
	cout << "Enter k1: ";
	cin >> k[0];
	cout << "Enter k2: ";
	cin >> k[1];
	cout << "Enter k3: ";
	cin >> k[2];
	cout << "Enter k4: ";
	cin >> k[3];
	cout << "Enter d: ";
	cin >> d;
	cout << "Enter N: ";
	cin >> N;
	cout << "Enter 0 - for an indirect sign, 1 - for the main: ";
	cin >> type;
	runEvolutional(etalon, k, d, N, type);
}
