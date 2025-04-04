
#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_linalg.h>
using namespace std;

default_random_engine generator;
normal_distribution<double> distribution (0, 0.1);

int num = 2;
int CROSSBOUND = 100;
double majorAxis = 0.75;
double minorAxis = 0.75;
double PI = 3.14159;
double DIST = 36;
double THRESHOLD = exp (-majorAxis*majorAxis);
double epsilon, INF;
double varX, varY;

vector<double> areas, lengths, angles, eps;
FILE *out = fopen ("sketch_plot.txt", "w");

class Point {
public:

    Point (){}
    Point(double x, double y);

    Point operator+ (Point& first);

    Point operator- (Point& first);

    Point operator+= (Point& first);

    Point operator-= (Point& first);

    Point operator-();

    double length();

    double getX() const;

    double getY() const;

    bool operator< (const Point& other) const {
        return other.x < this->x || (other.x == this->x && other.y < this->y);
    }

    double x,y;

};

Point drone_start_B, drone_start_A;
vector<Point> gaussianCenter;
vector<Point> gaussianVar;

class Line {

public:
    Line(double m, double b);

    static Line buildByPoints(Point& start, Point& end);

    static Line buildByPointAndAngle(Point& start, double angle);

    double getM() const;

    double getB() const;

private:
    double m, b;
};


class PointUtil {

public:
    static double orientation(Point& one, Point& two, Point& three);

    static Point vector(double angle, double length);

    static Point perpendicular(Point& one, Point& two, double length, int orientation);

    static const int CLOCKWISE = 1;
    static const int COUNTERCLOCKWISE = -1;
};


class LineSegment {

public:
    LineSegment(Point &start, Point &end);

    LineSegment(const Line &line, const Point &start, const Point &end);

    LineSegment(const LineSegment& copySegment);

    ~LineSegment();

    double length();

    Line getLine();

    Point getStart();

    Point* getStartPtr();

    Point getEnd();

    Point* getEndPtr();

    Line line;
    Point *start, *end;
};

typedef pair<Point, int> CrossData;
typedef pair<double, double> Pair;

class Ellipse {

public:
    Ellipse (){}
    Ellipse(Point centergiven, double x_radius, double y_radius)
    {
        center = centergiven;
        radius_x = x_radius;
        radius_y = y_radius;
        _size = M_PI * radius_x * radius_y;
    }

    bool inside(const Point &vector);

    double size();

    LineSegment segmentIntersections(LineSegment &segment);

    LineSegment intersections(Line &line);

    bool crosses(LineSegment &segment);

    bool crossesEdge(LineSegment &segment);

    Point getCross(LineSegment &segment) ;

    double edgeGradient(Point& point);

    Point getCenter() {
        return center;
    }

    double getXRadius() {
        return radius_x;
    }

    double getYRadius() {
        return radius_y;
    }
    
    Point center;
    double radius_x, radius_y;
    double _size;
};


/*
 FIX ME:
    CREATE A CRITICAL PATH CLASS CONTAINING:
    THE INITIAL CRITICAL POINT ON THE OUTER CONTOUR
    THE DIFFERENCE BETWEEN VALUES OF SUCCESSIVE CONTOURS
    THE STARTING AND ENDPOINT OF EACH SUCCESSIVE CONTOURS
 */

struct criticalPath{
    Point stcriticalPoint;
    double diffepsilon;
    vector<Pair> contourInfo;
    
    CrossData getCross (LineSegment segment, double nabla, double alpha, double dist, int inside)
    {
       Point init = *segment.start;
       Point fini = *segment.end;
       Point curr = init;
        
       Pair gradient[2];
       
       Point motion = PointUtil::vector (nabla + alpha, dist/100);
       
        //fix me
       // bool cross = checkCross ();
        
        //fix me: find crossing point using contourInfo.
       
      // cout <<"Exception did not cross!"<<endl;
      // exit (0);
       
       return CrossData ( Point (0,0), 0 );
   }
};

double get_dist (Point A, Point B)
{
    return sqrt ((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

double changeGradient (double angleA, double angleB)
{
    return abs (angleA - angleB); //check this again for bug.
}

//FIX ME

bool checkifSource(Point S){
    /*
     FIX ME: fill this in, is S a local maxima??
     */
    return false;
}

void reverse (double &gradient)
{
    if (gradient < 0)
        gradient = gradient + PI;
    else if (gradient > 0)
        gradient = gradient - PI;
    else
        gradient = 0 ;
    return ;
}

int checkin (Point P){
    double A = majorAxis;
    double B = minorAxis;
    return (((P.x * P.x) / (A * A) + (P.y * P.y) / (B*B))  <= 1);
}

vector<double> get (vector<double> a, vector<double> b)
{
    for (int i = 0;i < a.size(); i ++)
        a[i] -= b[i];
    return a;
}

vector<vector<double> > prod (vector<vector<double> > A, vector<vector<double> > B)
{
    vector<vector<double> > product;
    
    int m = A.size();
    int n = A[0].size();
    
    for (int i = 0;i < m; i ++)
    {
        vector<double> row;
        for (int j = 0; j < B[0].size(); j ++)
        {
            double val = 0;
            for (int k = 0 ; k < n; k++)
                val += A[i][k] * B[k][j];
            row.push_back (val);
        }
        product.push_back (row);
    }
    
    return product;
}

vector<double> prod (vector<vector<double> > A, vector<double> x)
{
    vector<double> product;
    
    for (int i = 0;i < A.size(); i ++)
    {
        double val = 0.0;
        for (int j = 0;j < A[i].size() ; j ++)
            val += A[i][j] * x[j];
        product.push_back (val);
    }
    
    return product;
}

vector<double> prod (vector<double> A, vector<double> x)
{
    vector<double> product;
    
    for (int i = 0;i < A.size(); i ++)
        A[i] = A[i] * x[0];
    
    product = A;
    
    return product;
}


double dot_prod (vector<double> A, vector<double> x)
{
    double total = 0;
    
    for (int i = 0;i < A.size(); i ++)
        total += A[i] * x[i];
        
    return total;
}

vector<vector<double> > transpose (vector<vector<double> > A)
{
    vector<vector<double> > matrix (A[0].size(), vector<double> (A.size(),0)) ;
    for (int i = 0;i < A.size(); i ++)
    {
        for (int j = 0;j < A[0].size(); j ++)
            matrix[j][i] = A[i][j];
    }
    return matrix;
}

vector<vector<double> > transpose (vector<double> x)
{
    vector<vector<double> > matrix ;
    matrix.push_back (x);
    return matrix;
}

double inner_product(vector<double> vec)
{
    return prod (transpose(vec), vec)[0];
}

vector<double> get_gradient (vector<vector<double> > A, vector<double> x, vector<double> b)
{
    vector<double> gradient;
    int m = A.size();
    int n = x.size();
    
    for (int i = 0 ; i < n; i++)
    {
        vector<double> row;
        for (int j  = 0; j < m; j ++)
            row.push_back (A[j][i]);
        vector<double> temp;
        temp = get (prod(A,x),b);
        double val = 2 * dot_prod (row, temp);
        gradient.push_back (val);
    }
    
    return gradient ;
}

vector<double> gradient_descent_convex(int id, vector<vector<double> > A, vector<double> b, vector<Point> points, vector<double> fval, vector<double> last_gradient)
{
    int iter = 0;
    vector<double> x, y, z;
    vector<vector<double> > X;
    
    z = vector<double> (A[0].size(), 0);

    vector<double> eta;
    eta.push_back (PI/1000);
    double val ;
  //  x.push_back (-2*(points[1].x-gaussianCenter[id].x) * fval[1]);
  //  x.push_back (-2*(points[1].y-gaussianCenter[id].y) * fval[1]);
   // double norm = sqrt (x[0] * x[0] + x[1]*x[1]);
   // x[0]/=norm; x[1]/=norm;
    
 //   return x;/*
    x = last_gradient;
  
    while (iter < 1000)
    {
        y = get (x, prod (get_gradient(A,x,b), eta));
        x = y;
        X.push_back (x);
        iter++;
    }
    
    for (int i = 0;i < X.size(); i ++)
        z = get (z,X[i]);
    
    vector<double> temp;
    temp.push_back (-1.0/X.size());
    z = prod (z, temp);
    
    vector<double> err;
    
    err = get(b, prod (A,z));
    
 //   cout<<"Error is "<<inner_product (err)<<endl;
    
    return z;
}

double gaussian (Point input, int gaussianId)
{
    double X = input.x - gaussianCenter[gaussianId].x;
    double Y = input.y - gaussianCenter[gaussianId].y;
    double varX = gaussianVar[gaussianId].x;
    double varY = gaussianVar[gaussianId].y;
    return exp ( - (X * X / (2*varX)) - (Y * Y/(2*varY)));
}

vector<double> getGaussian (vector<Point> points)
{
    vector<double> fval;
    
    for (int i = 0; i < points.size(); i ++)
    {
        double sum = 0;
        for (int j = 0; j < gaussianCenter.size(); j++)
            sum = sum + gaussian (points[i], j);
        fval.push_back (sum);
    }
    
    return fval;
}

double getGaussian (Point point)
{
    double fval;
    double sum = 0;
    for (int j = 0; j < gaussianCenter.size(); j++){
        sum = sum + gaussian (point, j);
    } 
    fval = sum;
    return fval;
}

vector<double> get_Gaussian_vector (vector<Point> points, int id)
{
    vector<double> x;
    vector<double> fval = getGaussian (points);
    //below is buggy, does not incorporate variance
    x.push_back (-2*(points[1].x-gaussianCenter[id].x) * fval[1]);
    x.push_back (-2*(points[1].y-gaussianCenter[id].y) * fval[1]);
    return x;
}


double gradient_modulo (double gradient)
{
    if (gradient > PI)
        gradient = gradient - 2*PI;
    else if (gradient < -PI)
        gradient = gradient + 2*PI;
    return gradient;
}

vector<vector<double> > inverse (vector<vector<double> > A)
{
    double a , b , c , d;
    a = A[0][0] ;
    b = A[0][1] ;
    c = A[1][0] ;
    d = A[1][1] ;
    
    double det = a*d - b*c;
    
    if (abs(det) < 1e-9)
    {
        cout <<"Exception non-invertible!"<<endl;
        cout << a << " " << b << endl;
        cout << c << " " << d << endl;
        exit(0);
    }
    
    A[0][0] = d/det;
    A[0][1] = -b/det;
    A[1][0] = -c/det;
    A[1][1] = a/det;
    
    return A ;
}

vector<double> solve (vector<double> equation)
{
    vector<double> solutions;
    double discriminant = 0;
    double a = equation[0];
    double b = equation[1];
    double c = equation[2];
    
    discriminant = b*b - 4*a*c;
    
    if (discriminant < 0)
        return solutions;
    
    solutions.push_back ( (sqrt (discriminant) - b) / (2*a));
    solutions.push_back ( (-sqrt (discriminant) - b) / (2*a));
    
    return solutions;
}

vector<double> polynomial_product(const vector<double>& p1, const vector<double>& p2) {
    int degree1 = p1.size();
    int degree2 = p2.size();
    vector<double> result(degree1 + degree2 - 1, 0);

    for (int i = 0; i < degree1; ++i) {
        for (int j = 0; j < degree2; ++j) {
            result[i + j] += p1[i] * p2[j];
        }
    }

    return result;
}

vector<double> getnullspace(vector<vector<double>> A) {
    int rows = A.size();
    int cols = A[0].size();

    if (rows != cols || (rows != 2 && rows != 3)) {
        cerr << "Only 2x2 or 3x3 matrices are supported!" << endl;
        exit(1);
    }

    vector<double> nullspace;

    if (rows == 2) {
        double a1 = A[0][0];
        double b1 = A[0][1];
        double a2 = A[1][0];
        double b2 = A[1][1];

        double x = 0, y = 0;
        if (abs(b1) > 1e-9) {
            y = -a1 / b1;
            x = 1;
        } else if (abs(b2) > 1e-9) {
            y = -a2 / b2;
            x = 1;
        }

        double norm = sqrt(x * x + y * y);
        x /= norm;
        y /= norm;

        nullspace.push_back(x);
        nullspace.push_back(y);
    } else if (rows == 3) {
        Eigen::Matrix<double, 3, 3> M;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                M(i, j) = A[i][j];
            }
        }

        Eigen::FullPivLU<Eigen::Matrix<double, 3, 3>> lu(M);
        Eigen::Matrix<double, 3, 1> kernel = lu.kernel().col(0);

        for (int i = 0; i < 3; ++i) {
            nullspace.push_back(kernel(i));
        }
    }

    return nullspace;
}

vector<vector<double> > get_pseudo_inv (vector<vector<double> > A)
{
    vector<vector<double> > pseudo_inv = A;

    A = prod (transpose(A),A);
    /*
     Singular Value Decomposition
     */
    
    vector<double> equation;
    std::vector<double> eigenvalue;
    if (A.size() == 2) {
        equation = polynomial_product({-1, A[0][0]}, {-1, A[1][1]});
        equation[2] -= (A[1][0] * A[0][1]);
        eigenvalue = solve(equation);
    } else if (A.size() == 3) {
        Eigen::Matrix<double, 3, 3> M;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                M(i, j) = A[i][j];
            }
        }
        Eigen::EigenSolver<Eigen::Matrix<double, 3, 3>> solver(M);
        Eigen::VectorXcd ev = solver.eigenvalues();
        // Convert complex eigenvalues to a vector of doubles
        for (int i = 0; i < ev.size(); ++i) {
            if (ev[i].imag() < 1e-9) {
            // If the eigenvalue is real, store the real part
            eigenvalue.push_back(ev[i].real());
            } else {
            // If the eigenvalue is complex, store the magnitude
            eigenvalue.push_back(std::abs(ev[i]));
            }
        }
    }
    else
    {   
        cout << "Exception! Not 2x2 or 3x3 matrix!"<<endl;
        exit(0);
    }

    
 //   cout<<A[0][0]<<" "<<A[0][1]<<endl;
 //   cout<<A[1][0]<<" "<<A[1][1]<<endl;
    vector<vector<double> > V, U, sigma;
    
 //   sort (eigenvalue.begin(), eigenvalue.end());
 //   reverse (eigenvalue.begin(), eigenvalue.end());
    
    for (int i = 0; i < eigenvalue.size(); i ++){
        double egval = eigenvalue[i];
        
  //      cout<<egval<<" egval "<<endl;
        
        vector<vector<double> > tonullspace = A;
        for(int j = 0; j < A.size(); j ++)
        {
            tonullspace[j][j] -= egval;
        }
        
        vector<double> eigenvector = getnullspace (tonullspace);
        V.push_back (eigenvector);
    //    cout<<V[i][0]<<" egvec "<<V[i][1]<<endl;
        vector<double> temp (A.size(),0);
        if (eigenvalue[i]>=0)
            temp[i] = sqrt (eigenvalue[i]);
        else if (abs(eigenvalue[i]) < 1e-9)
            temp[i] = 0;
        else
        {
            cout<<"Negative eigenvalue! "<<eigenvalue[i]<< endl;
            exit(0);
        }
        sigma.push_back (temp);
    }
    

    U = pseudo_inv;
    U = prod (U, V);
    U = prod (U, sigma);
    
  /*  for (int i = 0;i < 2; i ++){
        double norm = sqrt(U[i][0]*U[i][0] + U[i][1]*U[i][1]);
        if (norm > 0){
            U[i][0]/=norm; U[i][1]/=norm;
        }
    }*/
    
   // U = transpose (U);
   // cout<<sigma[0][0]<<" "<<sigma[1][1]<<endl;
    for (int i = 0; i < 2; i ++){
        if (abs(sigma[i][i])>0){
            sigma[i][i] = 1/sigma[i][i] + 1e-9;
        }
    }
        
    pseudo_inv = V;
    pseudo_inv = prod (pseudo_inv, sigma);
    pseudo_inv = prod (pseudo_inv, transpose(U));
    
    return pseudo_inv;
}

vector<double> gradient_matrix_solver (vector<vector<double> > A, vector<double> b)
{
    
    vector<vector<double> > inv;
  //  inv = get_pseudo_inv (prod (transpose(A),A));
    inv = get_pseudo_inv (A);
  /*  cout<<"A below"<<endl;
    cout<<A[0][0]<<" "<<A[0][1]<<endl;
    cout<<A[1][0]<<" "<<A[1][1]<<endl;
    cout<<"Pseudo inverse below. "<<endl;
    cout<<inv[0][0]<<" "<<inv[0][1]<<endl;
    cout<<inv[1][0]<<" "<<inv[1][1]<<endl;
   */
 //   inv = prod (inv, transpose(A));
    vector<double> gradient = prod (inv, b);
  //  cout<<gradient[0]<<" gr "<<gradient[1]<<endl;
    return gradient;
}

//TODO: Note assumptions.......
vector<double> gradient_LSQ (vector<Point> points)
{
  //  for (auto P : points)
    //    cout<<P.x<<" P "<<P.y<<endl;
    vector<double> fval = getGaussian (points);
    vector<vector<double> > A;
    vector<double> b;
    
    b.push_back (fval[0] - fval[1]);
    b.push_back (fval[2] - fval[1]);
    
    vector<double> row;
    row.push_back (points[0].x - points[1].x);
    row.push_back (points[0].y - points[1].y);
    
    double gr[2]={0,1};
    
    if (abs(row[0]) < 1e-9)
        ;//   return get_Gaussian_vector(points,0);
    row[1] /= row[0] + 1e-9;
    b[0] /= row[0] + 1e-9;
    row[0] = 1;
    A.push_back (row);
//    cout << row[1] << " row " << b[0]<<endl;

    
    row.erase (row.begin(), row.end());
    row.push_back (points[2].x - points[1].x);
    row.push_back (points[2].y - points[1].y);
    
    if (abs(row[0]) < 1e-9)
        ;//    return get_Gaussian_vector(points,0);
    row[1] /= row[0] + 1e-9;
    b[1] /= row[0] + 1e-9;
    row[0] = 1;
    A.push_back (row);
    
 //   cout << row[1] << " row " << b[1]<<endl;
    return gradient_matrix_solver (A,b);
 //   return gradient_descent_convex (id, A,b, points, fval, last_gradient);

}


//FIX ME - COME UP WITH NEW NAME 
std::vector<double> concentration_gradient_LSQ(
    std::vector<Point> surroundingPoints,
    Point& pointToSolve) {

    // Ensure there are at least 3 points (including the point to solve)
    if (surroundingPoints.size() < 2) {
        std::cerr << "Insufficient points for gradient calculation." << std::endl;
        return {};
    }

    // Extract the coordinates and Gaussian values for the surrounding points
    double xB = pointToSolve.x;
    double yB = pointToSolve.y;
    double fb = getGaussian(pointToSolve);

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    for (auto& point : surroundingPoints) {
        double x = point.x;
        double y = point.y;
        double f = getGaussian(point);

        // Build the rows of matrix A and vector b
        A.push_back({x - xB, y - yB});
        b.push_back(f - fb);
    }
    // Compute the pseudo-inverse of A
    std::vector<std::vector<double>> A_pinv = get_pseudo_inv(A);

    // Solve for x = A_pinv * b
    std::vector<double> gradient = prod(A_pinv, b); 

    return gradient;
}

// Function to compute the Hessian matrix
std::vector<std::vector<double>> hessian_LSQ(
    std::vector<Point> surroundingPoints,
    Point& pointToSolve,
    std::vector<double>& gradientAtPoint) {

    double xE = pointToSolve.x;
    double yE = pointToSolve.y;
    double fxe = gradientAtPoint[0];
    double fye = gradientAtPoint[1];
    double fe = getGaussian(pointToSolve);

    if (surroundingPoints.size() < 3) {
        std::cerr << "Insufficient points for gradient calculation." << std::endl;
        return {};
    }

    int numPoints = surroundingPoints.size();
    std::vector<std::vector<double>> A(numPoints, std::vector<double>(3, 0.0));
    std::vector<double> b(numPoints, 0.0);

    for (size_t i = 0; i < surroundingPoints.size(); ++i) {
        double x = surroundingPoints[i].x;
        double y = surroundingPoints[i].y;
        double f = getGaussian(surroundingPoints[i]); 

        A[i][0] = std::pow(x - xE, 2) + 1e-9; //adding to prevent divide by zero
        A[i][1] = (x - xE) * (y - yE) + 1e-9;
        A[i][2] = std::pow(y - yE, 2) + 1e-9;
        b[i] = f - fe - fxe * (x - xE) - fye * (y - yE);
    }
    // Compute the pseudo-inverse of A
    std::vector<std::vector<double>> A_pinv = get_pseudo_inv(A); 

    // cout << "solving for x below "<<endl;

    // Solve for x = A_pinv * b
    if (A_pinv[0].size() != b.size()) {
        cerr << "Matrix dimensions do not match for multiplication!" << endl;
        cerr << "A_pinv size: " << A_pinv[0].size() << ", b size: " << b.size() << endl;
        exit(1);
    }
    std::vector<double> x = prod(A_pinv, b);

    // Construct the Hessian matrix
    std::vector<std::vector<double>> H(2, std::vector<double>(2, 0.0));
    H[0][0] = x[0];
    H[0][1] = x[1];
    H[1][0] = x[1];
    H[1][1] = x[2];

    return H;
}

// Function to calculate curvature
double calc_curvature(std::vector<double>& g, 
    std::vector<std::vector<double>>& H) {
    // Gradient components
    double Fx = g[0];
    double Fy = g[1];

    // Hessian components
    double Fxx = H[0][0];
    double Fxy = H[0][1];
    double Fyy = H[1][1];

    // Numerator of the curvature formula
    double num = -Fy * Fy * Fxx + 2 * Fx * Fy * Fxy - Fx * Fx * Fyy;

    // Denominator of the curvature formula
    double denom = std::pow(Fx * Fx + Fy * Fy, 1.5);

    // Handle division by zero
    if (std::abs(denom) < 1e-9) {
        return NAN;
    }

    // Return the curvature
    return num / denom;
}

// The calc_curvature_LSQ function is a wrapper that calculates 
// the curvature at a point using LSQ.
std::pair<std::vector<double>, double> calc_curvature_LSQ(
    std::vector<Point> surroundingPoints,
    Point& pointToSolve)
{
    // Calculate the gradient at pointB using three points.
    std::vector<double> gradB = concentration_gradient_LSQ(surroundingPoints, pointToSolve); 
 
    // Calculate the Hessian using nine points, the computed gradient, and the concentration function.
    std::vector<std::vector<double>> H = hessian_LSQ(surroundingPoints, 
                                                     pointToSolve, gradB);
 
    // Calculate the curvature based on the gradient and Hessian.
    double curvature = calc_curvature(gradB, H);
 
    // Return the computed gradient and curvature.
    return std::make_pair(gradB, curvature);
}

// Helper function to compute numerical gradient
std::vector<double> compute_gradient(const std::vector<double>& values) {
    size_t n = values.size();
    std::vector<double> gradient(n, 0.0);

    if (n < 2) {
        std::cerr << "Insufficient data points for gradient calculation." << std::endl;
        return gradient;
    }

    // Central difference for interior points
    for (size_t i = 1; i < n - 1; ++i) {
        gradient[i] = (values[i + 1] - values[i - 1]) / 2.0;
    }

    // Forward difference for the first point
    gradient[0] = values[1] - values[0];

    // Backward difference for the last point
    gradient[n - 1] = values[n - 1] - values[n - 2];

    return gradient;
}

// Function to compute the difference between consecutive elements
std::vector<double> compute_diff(const std::vector<double>& values) {
    size_t n = values.size();
    std::vector<double> diff(n - 1, 0.0);

    for (size_t i = 0; i < n - 1; ++i) {
        diff[i] = values[i + 1] - values[i];
    }

    return diff;
}

// Function to compute curvature derivatives
std::vector<double> compute_curvature_derivative(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 const std::vector<double>& k) {
    // Compute gradients
    std::vector<double> dx = compute_gradient(x);
    std::vector<double> dy = compute_gradient(y);
    std::vector<double> dk = compute_gradient(k);

    // Compute ds (arc length differences)
    std::vector<double> diff_x = compute_diff(x);
    std::vector<double> diff_y = compute_diff(y);
    size_t n = diff_x.size();
    std::vector<double> ds(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        ds[i] = std::sqrt(diff_x[i] * diff_x[i] + diff_y[i] * diff_y[i]);
    }

    // Compute ds_mid (midpoints of ds)
    std::vector<double> ds_mid(ds.size() + 1, 0.0);
    ds_mid[0] = ds[0];
    for (size_t i = 1; i < ds.size(); ++i) {
        ds_mid[i] = (ds[i - 1] + ds[i]) / 2.0;
    }
    ds_mid[ds.size()] = ds[ds.size() - 1];

    // Compute dk_ds
    std::vector<double> dk_ds(k.size(), 0.0);
    for (size_t i = 0; i < k.size(); ++i) {
        dk_ds[i] = std::pow(dk[i] / ds_mid[i], 2);
    }

    return dk_ds;
}

double getAngle (vector<double> A)
{
    double normA = sqrt (A[0] * A[0] + A[1] * A[1]);
    double angle = acos (abs(A[0])/normA);
    if (A[0] >= 0 && A[1] >= 0)
        return angle;
    if (A[0] < 0 && A[1] >= 0)
        return PI - angle;
    if (A[0] < 0 && A[1] < 0)
        return PI + angle;
    return -angle;
}

struct PLUME{
    vector<Ellipse> ovals;
    
    PLUME ()
    {
    }
    
    int pointinEllipse (Point P, Ellipse oval)
    {
        double A = majorAxis;
        double B = minorAxis;
        P.x -= oval.center.x;
        P.y -= oval.center.y;
        double leftHandside = (P.x * P.x) / (A * A) + (P.y * P.y) / (B*B);
        if (leftHandside <= 1)
            return 1;
        else
            return 0;
    }
    
    double concentration (Point P)
    {
        double sum = 0;
        
        for (int i = 0;i < num; i++){
            sum = sum + gaussian (P, i);
        }
        
        return sum ;
    }
    
    /*
     
    FIX ME:
     HERE WE NEED TO ADAPT THIS CROSSING FUNCTION FOR THE CASE OF CRITICAL PATH NOT ELLIPSE.
     
    CrossData getCross (LineSegment segment, double nabla, double alpha, double dist, int inside)
    {
        Point init = *segment.start;
        Point fini = *segment.end;
        Point curr = init;
        
        Point motion = PointUtil::vector (nabla + alpha, dist/100);
        
        for (int i = 0;i < 100; i ++){
            curr = curr + motion;
            if (inside && concentration (curr) < THRESHOLD){
                return CrossData (curr, 1);
            }
            else if (!inside && concentration (curr) > THRESHOLD)
                return CrossData (curr, 1);
        }
        
       // cout <<"Exception did not cross!"<<endl;
       // exit (0);
        
        return CrossData ( Point (0,0), 0 );
    }*/
  
    CrossData getCross (LineSegment segment, double nabla, double alpha, double dist, int inside)
    {
        Point init = *segment.start;
        Point fini = *segment.end;
        Point curr = init;
        
        Point motion = PointUtil::vector (nabla + alpha, dist/100);
        
        for (int i = 0;i < 100; i ++){
            curr = curr + motion;
            if (inside && concentration (curr) < THRESHOLD){
                return CrossData (curr, 1);
            }
            else if (!inside && concentration (curr) > THRESHOLD)
                return CrossData (curr, 1);
        }
        
       // cout <<"Exception did not cross!"<<endl;
       // exit (0);
        
        return CrossData ( Point (0,0), 0 );
    }
   /*
    CrossData getCross (LineSegment segment)
    {
        vector<int> candidates;
        for (int i = 0;i < num;  i++){
            int cross = pointinEllipse (*segment.start, ovals[i]) ^ pointinEllipse (*segment.end, ovals[i]);
            if (cross == 1)
            {
                candidates.push_back(i);
            //    cout<<"ovals[i] "<<candidates.back().x << " " << candidates.back().y<<endl;
            }
        }
        
     //   cout << candidates.size() << endl;
        
        for (int i = 0;i < candidates.size(); i ++)
        {
            int cnt = 1;
            for (int j = 0;j < num; j ++)
            {
                if (j != candidates[i])
                    cnt += pointinEllipse (ovals[candidates[i]].getCross(segment), ovals[j]);
            }
        //    cout << candidates.size()<<" tst "<< cnt << endl;
            if (cnt == 1)
                return CrossData (ovals[candidates[i]].getCross (segment), candidates[i]);
        }
        
     //   cout <<"Exception at crossing! "<<endl;
        return CrossData ( Point (0,0), 0 );
    }
    double edgeGradient (CrossData P)
    {
        return ovals[P.second].edgeGradient (P.first);
    }*/
};

int currentGaussian1;

class Drone {
    public :
    Point position;
    Point last;
    double nabla ; //just the gradient angle, not the slope
    int inside, numCross ;
    bool droneIn;
    double angleTurned;
    double distTraversed;
    int currentGaussian;
    vector<Point> polytope;
    
    Drone (){}
    Drone (Point P1, Point P2, int in, double nab, bool flag)
    {
        position = P1;
        last = P2;
        inside = in;
        nabla = nab;
        droneIn = flag;
        angleTurned = 0;
        distTraversed = 0;
        numCross = 0;
        currentGaussian = 0;
     //   polytope.push_back (position);
    }
    
    bool MoveDrone (double alpha, double dist, PLUME &plume, int callSource)
    {
        Point motion;
        Point nextPosition;
        vector<Point> points;
        points.push_back (position);
        
        motion = PointUtil::vector (nabla + alpha, dist);
        nextPosition = position + motion;
        
        points.push_back (nextPosition);
        LineSegment dronemotion = LineSegment (position, nextPosition);
        
        distTraversed += dist;
        
        if (abs(nextPosition.x) > 3 || abs(nextPosition.y) > 3)
            cout << "position exception! "<<endl;
        
       // bool cross = plume.crossesEdge (dronemotion);
        CrossData cd = plume.getCross (dronemotion, nabla, alpha, dist, inside);
      //  CrossData cd = plume.getCross (dronemotion);
        Point crossingPoint = cd.first;
        
        points.push_back (crossingPoint);
        
        swap (points[1], points[2]);
        int cross = cd.second; //changed to boolean if 1 then crossed, 0 did not cross.
     //   if (abs(crossingPoint.x) > 1e-9 || abs (crossingPoint.y) > 1e-9)
       //     cross = 1;
       // if (cross)
        //    cout<<cross<<  " cross "<<crossingPoint.x << " " <<crossingPoint.y<<endl;
        if (cross)
        {
            ++numCross;
            inside = inside ^ 1;
           // if (abs(alpha) > epsilon)
             //   cout<<"Testing alpha " << inside<<endl;
            vector<double> gradient_vector;
            
          //  if (cd.second != currentGaussian1)
            //    gradient_vector = get_Gaussian_vector( points, cd.second);
           // else
                gradient_vector = gradient_LSQ (points);
            
       //     last_gradient = gradient_vector;
       //     Point gradient_shift = PointUtil::vector (1.0/100,1);
       //     last_gradient[0] += gradient_shift.x;
       //     last_gradient[1] += gradient_shift.y;
            double angle = getAngle (gradient_vector);
            
           // if (!callSource)
             //   cout<<"Called from CrossPlume "<<endl;
            
            angle = gradient_modulo (angle);
            if (droneIn)
            {
          //      cout<<gradient_vector[0]<<" "<<gradient_vector[1] << " gradient vector,Points "<<crossingPoint.x<<" "<<crossingPoint.y<<endl;
             //   cout<<position.x<<" "<<position.y<<" "<<nextPosition.x<<" "<<nextPosition.y<<endl;
            }
            double gradient = angle + PI/2 ;
       //     cout<<angle<<" angle "<<points.size()<<" "<<gradient<< endl;
            gradient = gradient_modulo (gradient);
            
            Point checkPoint = PointUtil::vector (gradient, dist);
           
      //      cout<<"Checkpoint "<<checkPoint.x<<" "<<checkPoint.y<<endl;

            checkPoint = crossingPoint + checkPoint;
           
            
        //    cout<<"Checkpoint "<<checkPoint.x<<" "<<checkPoint.y<<endl;
            
            int orient ;
            if (inside)
                orient = PointUtil::CLOCKWISE;
            else
                orient = PointUtil::COUNTERCLOCKWISE;
           // else
             //   orient = (inside) ? (PointUtil::COUNTERCLOCKWISE) : (PointUtil::CLOCKWISE);
            
          //  if (callSource == 0)
        //        orient = orient ^ 1;
            
            
            if (PointUtil::orientation (position, crossingPoint, checkPoint) != orient)
                reverse (gradient);

            angleTurned += changeGradient (nabla + alpha, gradient);
            nabla = gradient;
            
        //    if (droneIn && !inside)
          //          cout << crossingPoint.x <<  " crossed here "<< crossingPoint.y <<  " " <<droneIn<<" inside: "<<inside<<" "<<" "<< gradient<< endl;
         }
        last = position;
        position = nextPosition;
        if (numCross)
            polytope.push_back (position);
      //  if (numCross == 1 && cross)
        //    polytope.push_back (cd.first);
    
        if (polytope.size() > 1)
        {
            int N = polytope.size();
       //     if (droneIn)
         //       fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", polytope[N-2].x, polytope[N-2].y, polytope[N-1].x, polytope[N-1].y);
            Point currtoinit = polytope.back() - polytope[0];
            
            /*
             FIX ME:
                ADD A TERMINATING CONDITION, WHENEVER EITHER DRONE OF THE PAIR
                REACHES A SOURCE, BUT HOW TO CHECK IF THEIR POSITION IS NEAR A SOURCE I.E. A LOCAL MAXIMA OF CONCENTRATION FUNCTION??
             */
            
            return (currtoinit.length() < INF) && (numCross > CROSSBOUND);
        }
        else
            return false ;
    }
};

// TODO
class DronePair {
    public:
        Drone droneA, droneB;
    
        DronePair() {}
        DronePair(Drone droneA, Drone droneB) : droneA(droneA), droneB(droneB) {}
   
        // FIXME
    
        bool MovePair(double alpha, double dist, PLUME &plume, int callSource) {
            bool moveA = droneA.MoveDrone (alpha, dist, plume, callSource);
            bool moveB = droneB.MoveDrone (alpha, dist, plume, callSource);
            return moveA & moveB;
        }
};

Point::Point(double x, double y) : x(x), y(y) {}

Point Point::operator+(Point &first) {
    return Point(x + first.x, y + first.y);
}

Point Point::operator+=(Point &first){
    this->x += first.x;
    this->y += first.y;
    return *this;
}

Point Point::operator-(Point &first) {
    return Point(x - first.x, y - first.y);
}

Point Point::operator-=(Point &first){
    this->x -= first.x;
    this->y -= first.y;
    return *this;
}

Point Point::operator-() {
    return Point(-this->x, -this->y);
}

double Point::length() {
    return sqrt(x * x + y * y);
}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

double PointUtil::orientation(Point& one, Point& two, Point& three) {
    double k=(two.getY() - one.getY())*(three.getX() - two.getX())-(two.getX() - one.getX()) * (three.getY() - two.getY());

    if(k>0) {
        return CLOCKWISE;
    } else {

        return COUNTERCLOCKWISE;
    }
}

Point PointUtil::vector(double angle, double length) {
    return Point(length * cos(angle), length * sin(angle));
}

Point PointUtil::perpendicular(Point &one, Point &two, double length, int orientation) {
    double delta_x = two.getX() - one.getX();
    double delta_y = two.getY() - one.getY();
    double angle = atan2(delta_y, delta_x);
    return vector(angle + (orientation * M_PI_2), length);
}


Line::Line(double m, double b) : m(m), b(b) {}

Line Line::buildByPoints(Point &start, Point &end) {
    double m = (end.getY() - start.getY()) / (end.getX() - start.getX() + 1e-9); //divide by zero case solved by 1e-9
    double b = start.getY() - (m * start.getX());

    return Line(m, b);
}

Line Line::buildByPointAndAngle(Point &start, double angle) {
    double m = tan(angle);
    double b = start.getY() - (m * start.getX());

    return Line(m, b);
}

double Line::getM() const {
    return m;
}

double Line::getB() const {
    return b;
}

LineSegment::LineSegment(Point &start, Point &end) : line(Line::buildByPoints(start, end)), start(new Point(start.getX(), start.getY())), end(new Point(end.getX(), end.getY())) {}

LineSegment::LineSegment(const Line &line, const Point &start, const Point &end) : line(line), start(new Point(start.getX(), start.getY())), end(new Point(end.getX(), end.getY())) {}

LineSegment::LineSegment(const LineSegment &copySegment): line(copySegment.line), start(new Point(copySegment.start->getX(), copySegment.start->getY())), end(new Point(copySegment.end->getX(), copySegment.end->getY())) {}

double LineSegment::length() {
    Point vector = (*end - *start);
    return vector.length();
}

Line LineSegment::getLine() {
    return line;
}

Point LineSegment::getStart() {
    return *start;
}

Point* LineSegment::getStartPtr() {
    return start;
}

Point LineSegment::getEnd() {
    return *end;
}

Point* LineSegment::getEndPtr() {
    return end;
}

LineSegment::~LineSegment() {
    delete start;
    delete end;
}

bool Ellipse::inside(const Point &vector) {
    return pow(((vector.getX() - center.getX()) / radius_x), 2) +
           pow(((vector.getY() - center.getY()) / radius_y), 2) <= 1;
}

double Ellipse::size() {
    return _size;
}

LineSegment Ellipse::segmentIntersections(LineSegment &segment) {
    Line line = segment.getLine();
    LineSegment intersectionSegment = intersections(line);

    Point start = intersectionSegment.getStart();
    Point end = intersectionSegment.getEnd();

    if(start.getX() > end.getX()) {
        Point tmp = start;
        start = end;
        end = tmp;
    }

    Point segmentStart = segment.getStart();
    Point segmentEnd = segment.getEnd();

    if(segmentStart.getX() > segmentEnd.getX()) {
        Point tmp = segmentStart;
        segmentStart = segmentEnd;
        segmentEnd = tmp;
    }

    if(segmentStart.getX() > start.getX()) {
        start = segmentStart;
    }
    if(segmentEnd.getX() < end.getX()) {
        end = segmentEnd;
    }
    if(start.getX() > end.getX()) {
        start = end;
    }

    return LineSegment(segment.getLine(), start, end);
}

LineSegment Ellipse::intersections(Line &line) {
    double rx = radius_x * radius_x;
    double ry = radius_y * radius_y;

    double a = (1 / rx) + (line.getM() * line.getM() / ry);
    double b = (2 * line.getB() * line.getM() / ry) - (2 * center.getX() / rx) - (2 * center.getY() * line.getM() / ry);
    double c = (line.getB() * line.getB() / ry) - (2 * line.getB() * center.getY() / ry) + (center.getX() * center.getX() / rx) + (center.getY() * center.getY() / ry) - 1;

    // Solution using Quadratic equation -b +- sqrt(b^2 - 4ac)/2a
    // where ax^2 + bx + c = 0
    double discriminant = pow(b,2) - (4 * a * c);

    if (discriminant > 0){
        double x1 = ((-b) + sqrt(discriminant)) / (2 * a);
        double x2 = ((-b) - sqrt(discriminant)) / (2 * a);

        double y1 = (line.getM() * x1) + line.getB();
        double y2 = (line.getM() * x2) + line.getB();

        return LineSegment(line, Point(x1, y1), Point(x2, y2));
    } else {
        return LineSegment(line, Point(0, 0), Point(0, 0));
    }
}

bool Ellipse::crosses(LineSegment &segment) {
    Line line = segment.getLine();

    double rx = radius_x * radius_x;
    double ry = radius_y * radius_y;

    double a = (1 / rx) + (line.getM() * line.getM() / ry);
    double b = (2 * line.getB() * line.getM() / ry) - (2 * center.getX() / rx) - (2 * center.getY() * line.getM() / ry);
    double c = (line.getB() * line.getB() / ry) - (2 * line.getB() * center.getY() / ry) + (center.getX() * center.getX() / rx) + (center.getY() * center.getY() / ry) - 1;

    // Solution using Quadratic equation -b +- sqrt(b^2 - 4ac)/2a
    // where ax^2 + bx + c = 0
    double discriminant = pow(b,2) - (4 * a * c);

    if (discriminant > 0){
        double x1 = ((-b) + sqrt(discriminant)) / (2 * a);
        double x2 = ((-b) - sqrt(discriminant)) / (2 * a);

        return (segment.getStart().getX() < x1 && segment.getStart().getX() > x2)
            || (segment.getEnd().getX() < x1 && segment.getEnd().getX() > x2);
    } else {
        return false;
    }
}


bool Ellipse::crossesEdge(LineSegment &segment) {
    Line line = segment.getLine();
    LineSegment intersectionSegment = intersections(line);

    Point start = intersectionSegment.getStart();
    Point end = intersectionSegment.getEnd();

    Point segmentStart = segment.getStart();
    Point segmentEnd = segment.getEnd();
    
    //check if discriminant is negative
    if (start.x == 0 && start.y == 0 && end.x == 0 && end.y == 0)
        return false;

    if(segmentStart.getX() > segmentEnd.getX()) {
        Point tmp = segmentStart;
        segmentStart = segmentEnd;
        segmentEnd = tmp;
    }

    return (start.getX() > segmentStart.getX() && start.getX() < segmentEnd.getX()) ||
            (end.getX() > segmentStart.getX() && end.getX() < segmentEnd.getX());
}


Point Ellipse::getCross(LineSegment &segment) {
    Line line = segment.getLine();
    LineSegment intersectionSegment = intersections(line);

    Point start = intersectionSegment.getStart();
    Point end = intersectionSegment.getEnd();

    Point segmentStart = segment.getStart();
    Point segmentEnd = segment.getEnd();
    
    
    //check if discriminant is negative
    if (abs(start.x) <= 1e-9 && abs(start.y) <= 1e-9 && abs(end.x) <= 1e-9 && abs(end.y) <= 1e-5){
        cout<<"Exception with discriminant!!"<<endl;
        return Point (0,0);
    }
    
    if (abs(get_dist (start, segmentStart) + get_dist (start, segmentEnd) - get_dist (segmentStart,segmentEnd)) <= 1e-5)
        return start;

    
    if (abs(get_dist (end, segmentStart) + get_dist (end, segmentEnd) - get_dist (segmentStart,segmentEnd)) <= 1e-5)
        return end;

    cout<<"Exception! at Ovals "<<segmentStart.x << " "<<segmentStart.y<<" " <<segmentEnd.x << " " <<segmentEnd.y<< endl;
    cout<<"Ovals cont. "<<start.x<<" "<<start.y<<" "<<end.x<<" "<<end.y << endl;
    return Point(0, 0);
}

double Ellipse::edgeGradient(Point& point) {
    return atan2(-pow(radius_y, 2) * (point.getX()-center.getX()), (pow(radius_x, 2) * (point.getY()-center.getY())));
}

void Sync (Drone &A, Drone &B, double alpha, double dist, PLUME &plume)
{
    Point v;
    
    if (alpha > 0)
        v = PointUtil::vector (A.nabla + alpha - PI/2, DIST * epsilon);
    else
        v = PointUtil::vector (A.nabla + alpha + PI/2, DIST * epsilon);
    
    Point nextPosition = A.position + v;
    B.last = nextPosition;
    B.position = nextPosition;
    LineSegment dronemotion = LineSegment (B.last, B.position);
 //   if (abs(B.nabla - A.nabla) < 1e-9 &&  plume.getCross(dronemotion, B.nabla, alpha, dist, B.inside).second == 1)
   //     B.inside = B.inside ^ 1;
    B.nabla = A.nabla;
   // B.inside = A.inside;
    B.polytope.push_back (nextPosition);
    //ignoring angle turned during sync.
    
    return ;
}

string check (Point P){
    double A = majorAxis;
    double B = minorAxis;
    return (((P.x * P.x) / (A * A) + (P.y * P.y) / (B*B))  <= 1) ? "Inside " : "Outside ";
}


void print_data(Drone A, Drone B)
{
    int numPoints = A.polytope.size();
    
    
    fprintf (out, "Pen b\n");
    
    for (int i = 0;i < numPoints; i ++)
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", A.polytope[i].x, A.polytope[i].y, A.polytope[(i+1)%numPoints].x, A.polytope[(i+1)%numPoints].y);
  
    numPoints = B.polytope.size();
  
    cout <<numPoints << endl;
  
//    fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[0].center.x, plume.ovals[0].center.y, R, R/2);
 //   fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[1].center.x, plume.ovals[1].center.y, R, R/2);

    fprintf (out, "Pen r\n");
    
    for (int i = 0;i < numPoints; i ++)
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", B.polytope[i].x, B.polytope[i].y, B.polytope[(i+1)%numPoints].x, B.polytope[(i+1)%numPoints].y);
  
    return ;
}

bool CrossPlume (Drone &A, Drone &B, double alpha, PLUME &plume)
{
    Point start_pos = A.last;
   // if (abs(start_pos.x) < 2 && abs(start_pos.y) < 2)
     //   cout <<"start pos ... "<< start_pos.x << " " << start_pos.y <<" "<<check(A.position)<< endl;

    int crossing = A.inside;
    bool orient = true ;
    double alphainitial = alpha;
    bool endHere = false;
    
    int iterate = 1000;
    
    do{
        endHere = endHere || A.MoveDrone (alpha, epsilon * epsilon, plume, 0);
     //   B.MoveDrone (alpha, epsilon * epsilon, plume, 0);
        
        if (PointUtil::orientation(A.last, A.position, start_pos) == PointUtil::CLOCKWISE && alpha > 0)
            orient = false ;
        if (PointUtil::orientation(A.last, A.position, start_pos) == PointUtil::COUNTERCLOCKWISE && alpha < 0)
            orient = false ;
      //  if (crossing == A.inside)
        //    Sync (A,B,alpha, epsilon * epsilon, plume);
      //  else
        //    Sync (A,B,alphainitial, epsilon * epsilon, plume);
        
        if (alpha > 0)
            alpha += epsilon;
        else
            alpha -= epsilon;
        A.angleTurned += abs (alphainitial);
     //   cout << "testing cross plume ... "<< A.position.x << " " << A.position.y << endl;
     //   cout << "testing cross plume ... "<< B.position.x << " " << B.position.y << endl;
        iterate--;
        if (iterate < 0){
            print_data (A,B);
            exit (0);
        }
       // cout<<crossing<<" "<<A.inside<<" "<<A.droneIn<<endl;
    }while (crossing == A.inside && orient && !endHere);
        
    if (crossing == A.inside && !endHere){
        A.polytope.pop_back ();
        A.position = A.last;
        A.angleTurned -= abs (alphainitial);
        double dx = start_pos.getX () - A.position.getX();
        double dy = start_pos.getY () - A.position.getY();
        double gradient = atan2 (dy, dx);
       
    //    A.position = A.last;
     //   B.position = B.last;
        
        Point d1 = start_pos - A.position;
        Point motion;
        if (d1.length() > epsilon*epsilon)
            motion = PointUtil::vector (gradient, epsilon * epsilon);
        else
            motion = PointUtil::vector (gradient, d1.length());
        motion = A.position + motion;
        Point d2 = start_pos - motion;
     //   cout<<A.position.x<<" " <<A.position.y<<" "<<d1.length()<< " distances " << d2.length()<<endl;

        if (d1.length() < d2.length()){
            reverse (gradient);
        }
        A.angleTurned += changeGradient (A.nabla + alpha, gradient);
        A.nabla = gradient;
     //   B.nabla = gradient;
     //   cout << start_pos.x << " starting here " << start_pos.y << " "<<gradient<<" "<<check (start_pos) <<" "<<check(A.position)<<" "<<A.inside<< endl;
     //   cout << "testing cross plume2 ... "<< A.position.x << " " << A.position.y << endl;
     //   cout << "testing cross plume2 ... "<< A.last.x << " " << A.last.y << endl;
        int iter = 0;
        while (crossing == A.inside && !endHere){
            d1 = start_pos - A.position;
            if (d1.length() > epsilon * epsilon)
            {
                endHere = endHere || A.MoveDrone (0, epsilon*epsilon, plume, 0);
           //     Sync (A,B,alphainitial, epsilon * epsilon, plume);

            //    B.MoveDrone (0, epsilon * epsilon, plume, 0);
             //   B.MoveDrone (0, epsilon*epsilon, plume);
            }
            else{
                endHere = endHere || A.MoveDrone (0, d1.length(), plume, 0);
            //    Sync (A,B,alphainitial, d1.length(), plume);

             //   B.MoveDrone (0, d1.length(), plume, 0);

            //    B.MoveDrone (0, d1.length()/2, plume);
            }
            
         //   if (abs(A.position.x) > 2 || abs(A.position.y) > 2)
           //     break ;
            ++iter;
            if (iter > 10000)
            {
                cout<<"Iterations exceeding ..."<<endl;
                print_data(A,B);
                exit (0);
            }
        //    Sync2 (A,B, alpha);
            //cout << "testing cross plume2 ... "<< A.position.x << " " << A.position.y << endl;
      //      cout << "testing cross plume2 ... "<< B.position.x << " " << B.position.y << endl;
        }
      //  B.nabla = A.nabla;
     //   Sync (A,B,alphainitial, epsilon * epsilon, plume);
    }
    
    return endHere;
}

double estimateArea (vector<Point> polygon)
{
    double area = 0 ;
    int n = polygon.size ();
    
    for (int i = 0;i < n; i ++)
        area += polygon[i].x * polygon [(i+1)%n].y - polygon[i].y* polygon[(i+1)%n].x;
    
    area /=2 ;
    
    return area ;
}


void sketch_algorithm (double alpha)
{
    cout << "Running Sketch Algorithm for Epsilon = " << epsilon << endl;
    
    vector<Ellipse> ell;
    
    for (int i = 0;i < num; i ++)
        ell.push_back (Ellipse(gaussianCenter[i], majorAxis, minorAxis));
    
    PLUME plume;
    
    /*
     FIX ME
     
     INITIALIZE THE ALGORITHM WITH DRONE PAIRS
     AND APPROPRIATE LOCATION OF THE PAIRS
    
     */
    
    plume.ovals = ell;
   // fprintf (out, "Ellipse (%lf,%lf) %lf %lf ", plume.ovals[0].center.x, plume.ovals[0].center.y, majorAxis, minorAxis);
    
    Drone B (drone_start_B, drone_start_B, 1, 0, false);
    Drone A (drone_start_A, drone_start_A, 1, 0, true);
    
    int TOKEN = 2;
    
    bool loopEnd = false ;
    
  //  fprintf (out, "Pen b\n");
    
    
    do{
        int iter = 0;
       //FIX ME: CHANGE A, B TO DRONE PAIRS AND IN THE FUNCTION CALL CROSS PLUME.
        while ((A.numCross + B.numCross == 0 || 1 == A.inside + B.inside) && !loopEnd)
        {
         //   cout << "testing TOKEN ... "<< A.position.x << " "<< A.position.y <<" 1"<<  endl;
            ++iter;
            if (iter > 10000)
            {
                cout<<"Iterations exceeding ..."<<endl;
                print_data(A,B);
                exit (0);
            }
            /*if (abs(A.position.x) > 2 || abs(A.position.y) > 2)
                break ;
            
            if (abs(B.position.x) > 2 || abs(B.position.y) > 2)
                break ;*/
            loopEnd = loopEnd || A.MoveDrone (alpha, epsilon, plume, 1);
      //      if (A.numCross > 1)
                loopEnd = loopEnd || B.MoveDrone (alpha, epsilon, plume, 1);
            //   Sync (A,B,alpha);
            
            
        }
        if (A.inside + B.inside != 1)
        {
            if (A.inside + B.inside == 0)
            {
                alpha = epsilon;
            //    B.nabla = A.nabla;

                loopEnd = loopEnd || CrossPlume (A,B, alpha, plume);
              //  if (A.inside + B.inside == 1)
                Sync (A,B,alpha, epsilon, plume);
            
                    B.nabla = A.nabla;
              
 //               cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
 //               cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
            }
            else
            {
            //    cout<<"Exception both inside!"<<endl;
          //      exit(0);
                alpha = -epsilon;
        //        Sync (B,A,alpha, epsilon, plume);
              //  A.nabla = B.nabla;

                loopEnd = loopEnd || CrossPlume (B,A, alpha, plume);
              //  if (A.inside + B.inside == 1)
                Sync (B,A,alpha, epsilon, plume);

                    A.nabla = B.nabla;
                
   //             cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
   //             cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
            }
        }
    }while (!loopEnd);
    
    cout <<"Initial crossing with A is " << A.polytope[0].x << " " << A.polytope[0].y<< endl;
    cout << "angle turned by A is " << " " << A.angleTurned << endl;
    cout << "distance traversed by A is "<< " " << A.distTraversed << endl;
    cout << "area estimated by A is " << " " << estimateArea (A.polytope) << endl;
    areas.push_back (estimateArea (A.polytope));
    lengths.push_back (A.distTraversed);
    angles.push_back (A.angleTurned);
    cout << "actual area is  " << PI * majorAxis * minorAxis << endl;
    
    print_data (A,B);
    
    return ;
}

void test_gen_sketch_additions(double alpha){
    cout << "Running Sketch Algorithm for Epsilon = " << epsilon << endl;
    
    vector<Ellipse> ell;
    
    for (int i = 0;i < num; i ++)
        ell.push_back (Ellipse(gaussianCenter[i], majorAxis, minorAxis));
    
    PLUME plume;
    
    plume.ovals = ell;
   // fprintf (out, "Ellipse (%lf,%lf) %lf %lf ", plume.ovals[0].center.x, plume.ovals[0].center.y, majorAxis, minorAxis);
    
    Drone B (drone_start_B, drone_start_B, 1, 0, false);
    Drone A (drone_start_A, drone_start_A, 1, 0, true);

    //create drone pair
    DronePair dronePair(A, B);
    // Generate eight points at equal angular intervals around the starting point
    vector<Point> points;
    Point startPoint = A.position; // Starting point of the drone
    // Generate points at 45-degree intervals (2 * PI / 8)
    for (int i = 0; i < 8; i++) {
        points.push_back(Point(startPoint.getX() + alpha*i, startPoint.getY() + alpha*i));
    }

    // test curvature function
    cout << "Testing curvature function..." << endl;
    std::pair<std::vector<double>, double> result = calc_curvature_LSQ(points, startPoint);
    cout << "Result of curvature function: " << result.first[0] << " " << result.first[1] << " " << result.second << endl;

}

void test_infrastructure()
{
    double alpha = 0.1;
    num = 2;
    CROSSBOUND = 100;
    majorAxis = 0.25;
    minorAxis = 0.25;
    DIST = sqrt (49);
    THRESHOLD = exp (-majorAxis*majorAxis);
    
    drone_start_B = Point (1,0);
    drone_start_A = Point (1,DIST * epsilon);
    
    gaussianVar.push_back (Point (5,5));
    gaussianVar.push_back (Point (1,5));
    
    gaussianCenter.push_back (Point (-2,0.5));
    gaussianCenter.push_back (Point (1,0));


    cout<<"Initializing test... "<<endl;
    cout<<"Epsilon "<<epsilon<<endl;
    cout<<"Initial direction "<<alpha<<endl;
    cout<<"Number of Gaussians "<<num<<endl;
    cout<<"Least difference between starting and end point "<<INF<<endl;
    cout<<"Minimum crossings before checking termination "<<CROSSBOUND<<endl;
   // cout<<"Major and Minor Axis of Gaussian "<<majorAxis<<" "<<minorAxis<<endl;
    cout<<"Minimum distance factor between drones and minimum distance "<<DIST<<" " <<DIST * epsilon<<endl;
    cout<<"Concentration THRESHOLD "<<THRESHOLD<<endl;
    cout<<"Drone B starting point "<<drone_start_B.x<<" "<<drone_start_B.y<<endl;
    cout<<"Drone B starting point "<<drone_start_A.x<<" "<<drone_start_A.y<<endl;
    
    cout<<"Centers of gaussians "<<endl;
    
    for (int i =0 ; i < num; i ++)
        cout<<gaussianCenter[i].x<<" "<<gaussianCenter[i].y<<endl;
        
    // sketch_algorithm(alpha);
    test_gen_sketch_additions(alpha);
}

int main()
{
    int i = 0;
    for (epsilon = 0.005; i < 20 ; epsilon += 0.001){
        INF = 4 * epsilon;
     //   if (epsilon == 0.007) continue;
        test_infrastructure ();
        eps.push_back (epsilon);
        ++i;
    }
    
    cout <<"areas ";
    for (int i = 0;i < areas.size(); i ++)
        cout << areas[i] << " ";
    
    cout << endl;
    
    
    cout <<"lengths ";
    for (int i = 0;i < lengths.size(); i ++)
        cout << lengths[i] << " ";
    
    cout << endl;
    
    
    cout <<"angles ";
    for (int i = 0;i < angles.size(); i ++)
        cout << angles[i] << " ";
    
    cout << endl;
    
    cout <<"epsilons ";
    for (int i = 0;i < eps.size(); i ++)
        cout << eps[i] << " ";
    
    cout << endl;
    
    return 0;
}


/*
 
 EPSILON = 0.05
 Point start_a (0.479583,0.0707107);
 Point start_b (0.5,0);
 angle turned by A is  97.4447
 distance traversed by A is  10.305
 
 EPSILON = 0.01
 Point start_a (0.499199,0.0141421);
 Point start_b (0.5,0);
 
 angle turned by A is  184.588
 distance traversed by A is  10.1943
 
 EPSILON = 0.005
 Point start_a (0.4998,0.00707107);
 Point start_b (0.5,0);
 angle turned by A is  63.9424
 distance traversed by A is  2.91295
 
 EPSILON = 0.001
 Point start_a (0.499992,0.00141421);
 Point start_b (0.5,0);
 
 angle turned by A is  185.67
 distance traversed by A is  2.929

 */
