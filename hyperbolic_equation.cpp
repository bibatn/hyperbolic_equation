#include "hyperbolic_equation.h"


double Functions::AnalyticalSolution(double x, double y, double z, double t) const
{
    return sin(M_PI * x / g.L_x) * sin(M_PI * y / g.L_y) * sin(2 * M_PI * z / g.L_z) * cos(a * t);
}

double Functions::Phi(double x, double y, double z) const
{
    return AnalyticalSolution(x, y, z, 0);
}






