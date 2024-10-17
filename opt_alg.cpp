#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try {
        double* p = new double[2]{0, 0};
        int i = 0;
        double x1 = x0 + d;

		double x[2] = { x0, x1 };
        
        if (ff(x[1], ud1, ud2) == ff(x[0], ud1, ud2)) {
            p[0] = x[0];
            p[1] = x[1];
            return p;
        }

		if (ff(x[1], ud1, ud2) > ff(x[0], ud1, ud2)) {
			d = -d;
			x[1] = x[0] + d;
			if (ff(x[1], ud1, ud2) >= ff(x[0], ud1, ud2)) {
				p[0] = x[0];
				p[1] = x[0] - d;
				return p;
			}
		}
        while (ff(x[1], ud1, ud2) <= ff(x[0], ud1, ud2)) {
        
            if (solution::f_calls > Nmax) {
                throw std::runtime_error("error");
            }

            i++;
            x[1] = x0 + pow(alpha, i) * d;
        }  

        if (d > 0) {
            p[0] = x[0];
            p[1] = x[1];
        } else {
            p[0] = x[1];
            p[1] = x[0];
        }
        return p;
    }
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		int k = 1;
    	std::vector<double> fi = {0, 1};

		while (fib[k] < (b - a) / epsilon) {
			k++;
			fi.push_back(fi[k - 1] + fi[k - 2]);
    	}

		double a[2] = {a, 0};
		double b[2] = {b, 0};
		double c[2] = {(b[0] - fi[k-1]) / (fi[k] * (b[0] - a[0])), 0};
		double d[2] = {a[0] + b[0] - c[0], 0};

		for (int i = 0; i < k - 3; i++) {
			if (ff(c[0], ud1, ud2) < ff(d[0], ud1, ud2)) {
				a[1] = a[0];
				b[1] = d[0];
			} else {
				b[1] = b[0];
				a[1] = c[0];
			}
			c[1] = b[1] - fi[k-i-2] / (fi[k-i-1] * (b[1] - a[1]));
			d[1] = a[1] + b[1] - c[1];
    	}
		Xopt.x = c[1]
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a0, double b0, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double c0 = (b0 + a0) / 2;
        double a[2] = {a0, 0}; 
        double b[2] = {b0, 0}; 
        double c[2] = {c0, 0}; 
        double d[2] = {0, 0};
		do {
			double l = ff(a[0], ud1, ud2) * (b[0] * b[0] - c[0] * c[0]) +
					   ff(b[0], ud1, ud2) * (c[0] * c[0] - a[0] * a[0]) + 
					   ff(c[0], ud1, ud2) * (a[0] * a[0] - b[0] * b[0]);
			double m = ff(a[0], ud1, ud2) * (b[0] - c[0]) + 
					   ff(b[0], ud1, ud2) * (c[0] - a[0]) + 
					   ff(c[0], ud1, ud2) * (a[0] - b[0]);
			if (m <= 0)
				throw string("error1");
			d[0] = d[1];
			d[1] = 0.5 * l / m;
			if(a[0] < d[1] && d[1] < c[0]) {
				if (ff(d[1], ud1, ud2) < ff(c[0], ud1, ud2)) {
					a[1] = a[0];
					c[1] = d[1];
					b[1] = c[0];
				} else {
					a[1] = d[1];
					c[1] = c[0];
					b[1] = b[1];
				}
			} else {
				if (c[0] < d[1] && d[1] < b[0]) {
					if (ff(d[1], ud1, ud2) < ff(c[0], ud1, ud2)) {
						a[1] = c[0];
						c[1] = d[1];
						b[1] = b[0];
					} else {
						a[1] = a[0];
						c[1] = c[0];
						b[1] = d[1];
					}
				} else {
					throw string("error2");
				}
				i++;
				if (solution::f_calls > Nmax) {
                	throw std::runtime_error("error3");
            	}
			}
		} while ( (b[0] - a[0]) < epsilon || abs(d[1] - d[0]) < gamma );

		Xopt.x = d[1];
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
