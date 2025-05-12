
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
#include <map>
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

    bool isNear(Point& other, double eps = 1e-6) {
        return std::abs(x - other.x) < eps && std::abs(y - other.y) < eps;
    }

    bool operator==(const Point& other) const {
        return std::abs(x - other.x) < 1e-9 && std::abs(y - other.y) < 1e-9;
    }

    double x,y;

};

Point drone_start_B, drone_start_A, drone_start_AB, drone_start_BB;
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
    LineSegment(Point start, Point end);

    LineSegment(const Line &line, const Point &start, const Point &end);

    LineSegment(const LineSegment& copySegment);

    // ~LineSegment();

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



double get_dist (Point A, Point B)
{
    return sqrt ((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}

double changeGradient (double angleA, double angleB)
{
    return abs (angleA - angleB); //check this again for bug.
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


vector<Point> removeDuplicatePoints(vector<Point>& points, double eps = 1e-9) {
    vector<Point> cleaned;
    if (points.empty()) return cleaned;

    cleaned.push_back(points[0]);

    for (size_t i = 1; i < points.size(); ++i) {
        if (!points[i].isNear(cleaned.back(), eps)) {
            cleaned.push_back(points[i]);
        }
    }
    return cleaned;
}

// FIX ME: use diff epsilon here?
vector<Point> getGaussianContours(double level,
                                  double stepSize,
                                  double xMin,
                                  double xMax,
                                  double yMin,
                                  double yMax) {

    vector<Point> levelSetPoints;

    for (double x = xMin; x < xMax; x += stepSize) {
        for (double y = yMin; y < yMax; y += stepSize) {
            // Grid corners
            Point p1(x, y);
            Point p2(x + stepSize, y);
            Point p3(x, y + stepSize);
            Point p4(x + stepSize, y + stepSize);

            double f1 = getGaussian(p1);
            double f2 = getGaussian(p2);
            double f3 = getGaussian(p3);
            double f4 = getGaussian(p4);

            // Check if the level set crosses the grid cell
            if ((f1 > level && f2 > level && f3 > level && f4 > level) ||
                (f1 < level && f2 < level && f3 < level && f4 < level)) {
                    // cout << "No crossing in this cell" << endl;
                    // cout << "cell: " << x << " " << y << endl;
                    // cout << f1 << " " << f2 << " " << f3 << " " << f4 << endl;
                continue; // No crossing
            }

            // Build the case index (bits: p1, p2, p4, p3)
            int idx = 0;
            if (f1 > level) idx |= 1;
            if (f2 > level) idx |= 2;
            if (f4 > level) idx |= 4;
            if (f3 > level) idx |= 8;

            // Interpolation lambda
            auto interp = [&](const Point& a, const Point& b, double fa, double fb) -> Point {
                double t = (level - fa) / (fb - fa);
                return Point(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
            };

            switch (idx) {
                case 0:
                case 15:
                    break;

                case 1:
                case 14:
                    levelSetPoints.push_back(interp(p1, p2, f1, f2));
                    levelSetPoints.push_back(interp(p1, p3, f1, f3));
                    break;

                case 2:
                case 13:
                    levelSetPoints.push_back(interp(p2, p1, f2, f1));
                    levelSetPoints.push_back(interp(p2, p4, f2, f4));
                    break;

                case 3:
                case 12:
                    levelSetPoints.push_back(interp(p1, p3, f1, f3));
                    levelSetPoints.push_back(interp(p2, p4, f2, f4));
                    break;

                case 4:
                case 11:
                    levelSetPoints.push_back(interp(p4, p2, f4, f2));
                    levelSetPoints.push_back(interp(p4, p3, f4, f3));
                    break;

                case 5: {
                    levelSetPoints.push_back(interp(p1, p2, f1, f2));
                    levelSetPoints.push_back(interp(p1, p3, f1, f3));
                    levelSetPoints.push_back(interp(p4, p2, f4, f2));
                    levelSetPoints.push_back(interp(p4, p3, f4, f3));
                    break;
                }

                case 6:
                case 9:
                    levelSetPoints.push_back(interp(p2, p1, f2, f1));
                    levelSetPoints.push_back(interp(p4, p3, f4, f3));
                    break;

                case 7:
                case 8:
                    levelSetPoints.push_back(interp(p1, p3, f1, f3));
                    levelSetPoints.push_back(interp(p4, p3, f4, f3));
                    break;

                case 10: {
                    levelSetPoints.push_back(interp(p1, p3, f1, f3));
                    levelSetPoints.push_back(interp(p2, p4, f2, f4));
                    levelSetPoints.push_back(interp(p1, p2, f1, f2));
                    levelSetPoints.push_back(interp(p3, p4, f3, f4));
                    break;
                }
            }
        }
    }




    vector<Point> cleaned = removeDuplicatePoints(levelSetPoints);

    if (cleaned.size() < 2) {
        cout << "few level set points: " << cleaned.size() << endl;
        return cleaned;
    }

    // cout << "Level set points: " << cleaned.size() << endl;
    // for (const auto& point : cleaned) {
    //     cout << point.x << " " << point.y << endl;
    // }

    return cleaned;
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
   
    // getting negatives discriminants that should be considered zero
    // if considered zero, we have one solution 
    if (abs(discriminant)  < 1e-9){ 
        solutions.push_back ((-b) / (2*a));
        solutions.push_back ((-b) / (2*a));
        return solutions;
    }else if (discriminant < 0){
        return solutions;
    }else{
        solutions.push_back ( (sqrt (discriminant) - b) / (2*a));
        solutions.push_back ( (-sqrt (discriminant) - b) / (2*a));
        
        return solutions;
    }
    
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

double determinant(const std::vector<std::vector<double>>& M) {
    int n = M.size();

    if (n == 2) {
        // 2x2 determinant: ad - bc
        return M[0][0] * M[1][1] - M[0][1] * M[1][0];
    } else if (n == 3) {
        // 3x3 determinant using cofactor expansion
        return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
             - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
             + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    } else {
        std::cerr << "Only 2x2 and 3x3 matrices are supported.\n";
        exit(1);
    }
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

        // double det = determinant(A);
        // if(abs(det) > 1e-9){
        //     continue;
        // }

        vector<double> eigenvector = getnullspace(tonullspace);

        // If nullspace fails (zero vector, all zeros, or NaNs), use a fallback
            bool invalid = eigenvector.empty() || std::all_of(eigenvector.begin(), eigenvector.end(), [](double v) {
                return v == 0.0 || std::isnan(v) || std::isinf(v);
            });

        if (invalid) {
            eigenvector = std::vector<double>(A.size(), 0.0);
            eigenvector[i] = 1.0; // fallback: unit vector in ith direction
        }

        // Normalize
        double norm = std::sqrt(std::inner_product(eigenvector.begin(), eigenvector.end(), eigenvector.begin(), 0.0));
        for (double& val : eigenvector) val /= norm;
    
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

struct CubicSpline {
    std::vector<double> x, a, b, c, d;  

    void setPoints(const std::vector<double>& xs,
                    const std::vector<double>& ys) {
        int n = xs.size();
        x = xs;
        a = ys;

        std::vector<double> h(n-1), alpha(n-1);
        for (int i = 0; i < n-1; ++i)
            h[i] = xs[i+1] - xs[i];

        for (int i = 1; i < n-1; ++i)
            alpha[i] = 3.0*(a[i+1] - a[i]) / h[i]
                     - 3.0*(a[i]   - a[i-1]) / h[i-1];

        c.assign(n, 0.0);
        std::vector<double> l(n), mu(n), z(n);
        l[0] = 1.0;  mu[0] = z[0] = 0.0;
        for (int i = 1; i < n-1; ++i) {
            l[i] = 2*(xs[i+1] - xs[i-1]) - h[i-1]*mu[i-1];
            mu[i] = h[i] / l[i];
            z[i]  = (alpha[i] - h[i-1]*z[i-1]) / l[i];
        }
        l[n-1] = 1.0;  z[n-1] = c[n-1] = 0.0;

        b.assign(n-1, 0.0);
        d.assign(n-1, 0.0);
        for (int j = n-2; j >= 0; --j) {
            c[j] = z[j] - mu[j]*c[j+1];
            double hi = h[j];
            b[j] = (a[j+1]-a[j])/hi - hi*(c[j+1]+2*c[j])/3.0;
            d[j] = (c[j+1] - c[j]) / (3.0*hi);
        }
    }

    double operator()(double xi) const {
        int i = std::upper_bound(x.begin(), x.end(), xi) - x.begin() - 1;
        if (i < 0)            i = 0;
        else if (i >= b.size()) i = b.size()-1;
        double dx = xi - x[i];
        return a[i]
             + b[i]*dx
             + c[i]*dx*dx
             + d[i]*dx*dx*dx;
    }

    double getDerivative(double xi) const {
        int i = std::upper_bound(x.begin(), x.end(), xi) - x.begin() - 1;
        if (i < 0)            i = 0;
        else if (i >= b.size()) i = b.size()-1;
        double dx = xi - x[i];
        return b[i]
             + 2.0*c[i]*dx
             + 3.0*d[i]*dx*dx;
    }
};


std::vector<double> curvature_gradient_LSQ(
    std::map<Point, double> surroundingCurvatures,
    std::pair<Point, double> crossingPointCurvature) {

    // Ensure there are at least 3 points (including the point to solve)
    if (surroundingCurvatures.size() < 2) {
        std::cerr << "Insufficient points for gradient calculation." << std::endl;
        return {};
    }

    // Extract the coordinates and Gaussian values for the surrounding points
    Point pointToSolve = crossingPointCurvature.first;
    double xB = pointToSolve.x;
    double yB = pointToSolve.y;
    double kb = crossingPointCurvature.second;

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    for (auto& pointInfo : surroundingCurvatures) {
        Point point = pointInfo.first;
        double x = point.x;
        double y = point.y;
        double k = pointInfo.second;

        // Build the rows of matrix A and vector b
        A.push_back({x - xB, y - yB});
        b.push_back(k - kb);
    }
    // Compute the pseudo-inverse of A
    std::vector<std::vector<double>> A_pinv = get_pseudo_inv(A);

    // Solve for x = A_pinv * b
    std::vector<double> gradient = prod(A_pinv, b); 

    return gradient;
}

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
    // cout << "concentration at pts: " << fb << endl;

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    for (auto& point : surroundingPoints) {
        double x = point.x;
        double y = point.y;
        double f = getGaussian(point);
        // cout << "concentration at surps " << f << endl;

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

        A[i][0] = 0.5 * std::pow(x - xE, 2); 
        A[i][1] = (x - xE) * (y - yE);
        A[i][2] = 0.5 * std::pow(y - yE, 2);
        b[i] = f - fe - fxe * (x - xE) - fye * (y - yE);
    }

    std::vector<std::vector<double>> A_pinv = get_pseudo_inv(A); 

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
    double Fx = g[0];
    double Fy = g[1];

    double Fxx = H[0][0];
    double Fxy = H[0][1];
    double Fyy = H[1][1];

    double num = -Fy * Fy * Fxx + 2 * Fx * Fy * Fxy - Fx * Fx * Fyy;

    double denom = std::pow(Fx * Fx + Fy * Fy, 1.5);

    return num / denom;
}

// The calc_curvature_LSQ function is a wrapper that calculates 
// the curvature at a point using LSQ.
double calc_curvature_LSQ(
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
    return curvature;
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
std::vector<double> compute_curvature_derivative(const std::vector<Point>& path,
                                                 const std::vector<double>& k) {

    // Extract x and y coordinates from the path
    size_t n = path.size();
    std::vector<double> x(n), y(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = path[i].x;
        y[i] = path[i].y;
    }

    // Compute gradients
    std::vector<double> dx = compute_gradient(x);
    std::vector<double> dy = compute_gradient(y);
    std::vector<double> dk = compute_gradient(k);
    
    // for (size_t i = 0; i < n; ++i) {
    //     cout << dx[i] << " dx " << dy[i] << " dy " << k[i] << " k " << endl;
    // }

    // Compute ds (arc length differences)
    std::vector<double> diff_x = compute_diff(x);
    std::vector<double> diff_y = compute_diff(y);
    n = diff_x.size();
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
        dk_ds[i] = dk[i] / ds_mid[i];
    }

    // Compute moving average of curvature derivatives
    size_t windowSize = 3; // Define the window size for the moving average
    std::vector<double> smoothedCurvatureDerivatives(k.size(), 0.0);

    for (size_t i = 0; i < k.size(); ++i) {
        double sum = 0.0;
        size_t count = 0;

        for (size_t j = i; j < std::min(i + windowSize, k.size()); ++j) {
            sum += dk_ds[j];
            ++count;
        }

        smoothedCurvatureDerivatives[i] = sum / count;
        // cout << smoothedCurvatureDerivatives[i] << " smoothedCurvatureDerivatives " << i << endl;
    }

    dk_ds = smoothedCurvatureDerivatives;

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

// class DronePair;
// class Drone;
// struct criticalPath;

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
    Point motion; // FIX ME ; can delete
    vector<Point> polytope;
    vector<double> lastContourGradient;
    vector<double> currentContourGradient;
    LineSegment currPath;
    vector<double> lastTangent;
    vector<double> currentTangent;
    int side;
    
    Drone() 
    : position(Point(0, 0)), 
      last(Point(0, 0)), 
      inside(0), 
      nabla(0.0), 
      droneIn(false), 
      angleTurned(0.0), 
      distTraversed(0.0), 
      numCross(0), 
      currentGaussian(0), 
      lastContourGradient(), 
      currentContourGradient(), 
      motion(Point(0,0)),
      currPath(Point(0,0), Point(0,0)),
      side(0) {}

    Drone(Point P1, Point P2, int side, double nab, bool flag) 
        : position(P1), 
        last(P2), 
        inside(0), 
        nabla(nab), 
        droneIn(flag), 
        angleTurned(0.0), 
        distTraversed(0.0), 
        numCross(0), 
        currentGaussian(0), 
        lastContourGradient(), 
        currentContourGradient(), 
        motion(Point(0,0)),
        currPath(Point(0,0), Point(0,0)),
        side(side) {}
    
    bool MoveDrone (double alpha, double dist, PLUME &plume, int callSource)
    {
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
            
            
            return (currtoinit.length() < INF) && (numCross > CROSSBOUND);
        }
        else
            return false ;
    }

    // Fixed implementation
    void MoveDrone (double alpha, double dist, 
                                    double diffepsilon, int callSource)
    {
        Point nextPosition;
    
        // Calculate motion vector and next position
        motion = PointUtil::vector (nabla + alpha, dist);
        nextPosition = position + motion;
        
        // points.push_back (nextPosition);
        currPath = LineSegment (position, nextPosition);
        
        // Update distance traversed
        distTraversed += dist;
        
        // Check for position exception
        if (abs(nextPosition.x) > 5 || abs(nextPosition.y) > 5)
            cout << "position exception! "<<endl;

        last = position;
        position = nextPosition;
        polytope.push_back (position);

        cout << "moved one drone now getting contour grad" << endl;
        // Calculate gradient using criticalPath
        vector<Point> surroundingPoints;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                if (i == 0 && j == 0) continue; 
                if (i != 0 && j != 0) { 
                    // FIX ME: do we use diff epsilon here?
                    Point offset = Point(i * diffepsilon, j * diffepsilon);
                    surroundingPoints.push_back(position + offset);
                }
            }
        }

        vector<double> gradient = concentration_gradient_LSQ(surroundingPoints, position);
        // cout << "gradient size : " << gradient.size() << endl;
        // for(int i = 0; i < gradient.size(); i++){
        //     cout << "gradient: [" << i << "] " << gradient[i] << endl;
        // }
        if(lastContourGradient.size() > 0){
            lastContourGradient = currentContourGradient;
        }else{
            vector<Point> surroundingPointsLast;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if (i == 0 && j == 0) continue; 
                    if (i != 0 && j != 0) { 
                        // FIX ME: do we use diff epsilon here?
                        Point offset = Point(i * diffepsilon, j * diffepsilon);
                        surroundingPointsLast.push_back(last + offset);
                    }
                }
            }

            lastContourGradient = concentration_gradient_LSQ(surroundingPointsLast, last); 
        }
        currentContourGradient = gradient;

        cout << "current point: " << position.x << " " << position.y << endl;
    }

            // FIX ME: add function to deal with crossing and orientation
            void LearnGradient(double alpha, double dist, Point crossingPoint, Drone &otherDrone, vector<double> gradient_vector) {
                cout << "Learning gradient for drone pair" << endl;
                
                cout << "gradient vector: " << gradient_vector[0] << " " << gradient_vector[1] << endl;
    
                double angle = getAngle (gradient_vector);
    
    
                angle = gradient_modulo (angle);
                cout << "angle: " << angle << endl;
    
                double gradient = angle + PI/2 ;
                gradient = gradient_modulo (gradient);
                Point checkPoint = PointUtil::vector (gradient, dist);
                checkPoint = crossingPoint + checkPoint;
                        
                // FIX ME: double check this
                int orient;
                if (inside)
                    orient = PointUtil::CLOCKWISE;
                else
                    orient = PointUtil::COUNTERCLOCKWISE;
                
                Point curr =  position;
                if (PointUtil::orientation (curr, crossingPoint, checkPoint) != orient)
                    reverse (gradient);
    
                angleTurned += changeGradient(nabla + alpha, gradient);
                nabla = gradient;
                   
            }

};

class DronePair {
    public:
        Drone droneA, droneB;
        bool foundSource = false;
        int inside = 0; // Default initialization
        double nabla = 0.0; // Default initialization
        double angleTurned = 0.0; // Default initialization
        Point position = Point(0, 0); // Default initialization
        Point last = Point(0, 0); // Default initialization
        int numCross = 0; // Default initialization
        double pairDist = 0.0; // Default initialization
        vector<double> currentContourGradient; // Default initialization
        vector<double> lastContourGradient; // Default initialization

        DronePair() {}
        DronePair(Drone droneA, Drone droneB, double pairDist) : droneA(droneA), droneB(droneB), pairDist(pairDist) {}

        void setNabla(double nabla) {
            this->nabla = nabla;
            droneA.nabla = nabla;
            droneB.nabla = nabla;
        }

        double getNabla() {
            return nabla;
        }

        void setNumCross(int numCross) {
            this->numCross = numCross;
            droneA.numCross = numCross;
            droneB.numCross = numCross;
        }

        int getNumCross() {
            return droneA.numCross;
        }

        void setAngleTurned(double angleTurned) {
            this->angleTurned = angleTurned;
            droneA.angleTurned = angleTurned;
            droneB.angleTurned = angleTurned;
        }

        double getAngleTurned() {
            return angleTurned;
        }

        void setPosition(Point position) {
            this->position = position;
            droneA.position = position;
            // FIX ME
            Point dist = Point(sqrt(pairDist), sqrt(pairDist));
            droneB.position = position - dist;
        }
        
        Point getPosition() {
            return droneA.position;
        }

        void setLast(Point last) {
            this->last = last;
            droneA.last = last;
            // FIX ME
            Point dist =  Point(sqrt(pairDist), sqrt(pairDist));
            droneB.last = last - dist;
        }

        Point getLast() {
                return droneA.last;
            }

        void setInside(int inside) {
            this->inside = inside;
            droneA.inside = inside;
            droneB.inside = inside;
        }

        int getInside() {
            return droneA.inside;
        }

        void setSide(int side) {
            this->inside = side;
            droneA.side = side;
            droneB.side = side;
        }

        int getSide() {
            return droneA.side;
        }

        vector<double> getCurrentTangent() {
            return droneA.currentTangent;
        }

        vector<double> getLastTangent() {
            return droneA.lastTangent;
        }
    
        vector<double> getCurrentContourGradient() {
            return droneA.currentContourGradient;
        }

        vector<double> getLastContourGradient() {
            return droneA.lastContourGradient;
        }
         


        void MovePair(double alpha, double dist, double diffepsilon, int callSource) {
            droneA.MoveDrone(alpha, dist, diffepsilon, callSource); 
            droneB.MoveDrone(alpha, dist, diffepsilon, callSource);
            cout << "current point: " << getPosition().x << " " << getPosition().y << endl;

        }

        // FIX ME: add function to deal with crossing and orientation
        void LearnGradient(double alpha, double dist, Point crossingPoint, DronePair &otherDrone, vector<double> gradient_vector) {
            cout << "Learning gradient for drone pair" << endl;
            //FIX ME: may need to set differently for each drone 
            // setNumCross(getNumCross() + 1); 
            // setInside(inside ^ 1);
              
            // //learn curvature at the four contour points we took
            // // gradient at
            // std::map<Point, double> surroundingCurvatures;

            // // point at last contour  for this drone
            // Point pointToSolve = getLast();
            // vector<Point> surroundingPointsLastThis;
            // for (int i = -1; i <= 1; i++) {
            //     for (int j = -1; j <= 1; j++) {
            //         if (i == 0 && j == 0) continue; 
            //         if (i != 0 && j != 0) { 
            //             // FIX ME: do we use diff epsilon here?
            //             Point offset = Point(i * diffepsilon, j * diffepsilon);
            //             surroundingPointsLastThis.push_back(pointToSolve + offset);
            //         }
            //     }
            // }

            // cout << "solving for curvature at last contour point for this drone" << endl;
            
            // vector<double> contourGradient = getLastContourGradient();
            // std::vector<std::vector<double>> H = hessian_LSQ(surroundingPointsLastThis, 
            //     pointToSolve,contourGradient);

            // double curvature_lastThis = calc_curvature(contourGradient, H);

            // surroundingCurvatures[pointToSolve] = curvature_lastThis;

            // cout << "solving for curvature at current contour point for this drone" << endl;
            // // point at current cnotour for this drone      
            // vector<Point> surroundingPointsCurrentThis;
            // pointToSolve = getPosition(); 
            // for (int i = -1; i <= 1; i++) {
            //     for (int j = -1; j <= 1; j++) {
            //         if (i == 0 && j == 0) continue; 
            //         if (i != 0 && j != 0) { 
            //             // FIX ME: do we use diff epsilon here?
            //             Point offset = Point(i * diffepsilon, j * diffepsilon);
            //             surroundingPointsCurrentThis.push_back(pointToSolve + offset);
            //         }
            //     }
            // } 

            // H.clear();
            // contourGradient = getCurrentContourGradient();
            // H = hessian_LSQ(surroundingPointsCurrentThis, 
            //     pointToSolve, contourGradient);


            // double curvature_currentThis = calc_curvature(contourGradient, H);

            // surroundingCurvatures[pointToSolve] = curvature_currentThis;

            // cout << "solving for curvature at last contour point for other drone" << endl;
            // //point at last contour for other drone
            // vector<Point> surroundingPointsLastOther;
            // for (int i = -1; i <= 1; i++) {
            //     for (int j = -1; j <= 1; j++) {
            //         if (i == 0 && j == 0) continue; 
            //         if (i != 0 && j != 0) { 
            //             // FIX ME: do we use diff epsilon here?
            //             Point offset = Point(i * diffepsilon, j * diffepsilon);
            //             surroundingPointsLastOther.push_back(otherDrone.getLast() + offset);
            //         }
            //     }
            // } 

            // pointToSolve = otherDrone.getLast();

            // H.clear();
            // contourGradient = otherDrone.getLastContourGradient();
            // H = hessian_LSQ(surroundingPointsLastOther, 
            //     pointToSolve, contourGradient);

            // double curvature_lastOther = calc_curvature(contourGradient, H);

            // surroundingCurvatures[otherDrone.getLast()] = curvature_lastOther;

            // cout << "solving for curvature at current contour point for other drone" << endl;
            // //point at current contour for other drone
            // vector<Point> surroundingPointsCurrentOther;
            // for (int i = -1; i <= 1; i++) {
            //     for (int j = -1; j <= 1; j++) {
            //         if (i == 0 && j == 0) continue; 
            //         if (i != 0 && j != 0) { 
            //             // FIX ME: do we use diff epsilon here?
            //             Point offset = Point(i * diffepsilon, j * diffepsilon);
            //             surroundingPointsCurrentOther.push_back(otherDrone.getPosition() + offset);
            //         }
            //     }
            // } 

            // H.clear();
            // pointToSolve = otherDrone.getPosition();
            // contourGradient = otherDrone.getCurrentContourGradient();
            // H = hessian_LSQ(surroundingPointsCurrentOther, pointToSolve, contourGradient);

            // double curvature_currentOther = calc_curvature(contourGradient, H);

            // surroundingCurvatures[otherDrone.getLast()] = curvature_currentOther;

            // // get curvature at crossing point
            // vector<Point> surroundingPoints;
            // for (int i = -1; i <= 1; i++) {
            //     for (int j = -1; j <= 1; j++) {
            //         if (i == 0 && j == 0) continue; 
            //         if (i != 0 && j != 0) { 
            //             // FIX ME: do we use diff epsilon here?
            //             Point offset = Point(i * diffepsilon, j * diffepsilon);
            //             surroundingPoints.push_back(crossingPoint + offset);
            //         }
            //     }
            // } 

            // double curvature_crossPoint = calc_curvature_LSQ(surroundingPoints, crossingPoint);
            // std::pair<Point, double> crossingPointCurvatures(crossingPoint, curvature_crossPoint);

            // // now get gradient of curvatures at crossing point
            // vector<double> gradient_vector = curvature_gradient_LSQ(surroundingCurvatures,
            //                                                     crossingPointCurvatures);
            
            // // square
            // for(double& dk : gradient_vector){
            //     dk = dk * dk;
            // }

            // FIX ME:
            // vector<double> gradient_vector = cp.getGradientAtPoint(crossingPoint);
            
            cout << "gradient vector: " << gradient_vector[0] << " " << gradient_vector[1] << endl;

            double angle = getAngle (gradient_vector);


            angle = gradient_modulo (angle);
            cout << "angle: " << angle << endl;

            double gradient = angle + PI/2 ;
            gradient = gradient_modulo (gradient);
            Point checkPoint = PointUtil::vector (gradient, dist);
            checkPoint = crossingPoint + checkPoint;
                    
            // FIX ME: double check this
            int orient;
            if (inside)
                orient = PointUtil::CLOCKWISE;
            else
                orient = PointUtil::COUNTERCLOCKWISE;
            
            Point curr = getPosition();
            if (PointUtil::orientation (curr, crossingPoint, checkPoint) != orient)
                reverse (gradient);

            angleTurned += changeGradient (getNabla() + alpha, gradient);
            setNabla(gradient);
               
        }
};

struct criticalPath{
    Point stcriticalPoint;
    double diffepsilon;
    vector<Point> criticalPathPoints;
    double contourRes; // FIX ME: need to determine contour res?
    vector<double> lastCPCurvatures;
    static const int LEFT = 1;
    static const int RIGHT = 2;
    
    criticalPath(Point start, double epsilon, double res)
    : stcriticalPoint(start), diffepsilon(epsilon), contourRes(res) {
        // Initialize the critical path by finding the first critical point
        double level = getGaussian(start);
        vector<Point> initialContour = getGaussianContours(level, 
                                                           contourRes, 
                                                           start.x - DIST * epsilon, 
                                                           start.x + DIST * epsilon, 
                                                           start.y - DIST * epsilon, 
                                                           start.y + DIST * epsilon);
        if (!initialContour.empty()) {
            Point firstCriticalPoint = getCriticalPoint(initialContour);
            criticalPathPoints.push_back(firstCriticalPoint);
        } else {
            std::cerr << "Failed to initialize critical path: no contour points found." << std::endl;
        }
    }


    vector<Point> getCriticalPathPoints() {
        std::sort(criticalPathPoints.begin(), criticalPathPoints.end()); 
        vector<Point> cleaned = removeDuplicatePoints(criticalPathPoints, diffepsilon);

        return cleaned;
    }

    std::pair<CubicSpline, CubicSpline> getCriticalPathSpline() {
        // just want to sort just in case
        vector<Point> criticalPathPoints = getCriticalPathPoints();

        vector<double> xValues;
        vector<double> yValues;
        vector<double> t;
        double totalLength = 0.0;
        t.push_back(0.0);
        xValues.push_back(criticalPathPoints[0].x);
        yValues.push_back(criticalPathPoints[0].y);
        for (size_t i = 1; i < criticalPathPoints.size(); ++i) {
            double dx = criticalPathPoints[i].x - criticalPathPoints[i-1].x;
            double dy = criticalPathPoints[i].y - criticalPathPoints[i-1].y;
            totalLength += std::hypot(dx, dy);
            t.push_back(totalLength);
            xValues.push_back(criticalPathPoints[i].x);
            yValues.push_back(criticalPathPoints[i].y);
        }
        
        CubicSpline cp_spline_x;
        cp_spline_x.setPoints(t, xValues);

        CubicSpline cp_spline_y;
        cp_spline_y.setPoints(t, yValues);


        return std::make_pair(cp_spline_x, cp_spline_y);
    }

    vector<double> getGradientAtPoint(Point P){

        std::pair<CubicSpline, CubicSpline> cp_spline = getCriticalPathSpline();
        CubicSpline cp_spline_x = cp_spline.first;
        CubicSpline cp_spline_y = cp_spline.second;
    

        double tMin = cp_spline_x.x.front();
        double tMax = cp_spline_x.x.back();
        double bestT = tMin;
        double minDist = std::numeric_limits<double>::max();
    
        for (double t = tMin; t <= tMax; t += contourRes) { 
            double xt = cp_spline_x(t);
            double yt = cp_spline_y(t);
            double dist = std::hypot(xt - P.x, yt - P.y);
            if (dist < minDist) {
                minDist = dist;
                bestT = t;
            }
        }
    
        double dx = cp_spline_x.getDerivative(bestT);
        double dy = cp_spline_y.getDerivative(bestT);

        return {dx, dy};


    }

  

    vector<double> getContourTangent(Point A, Point B){
        Point current = A;
        vector<Point> surroundingPoints;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue; 
            if (i != 0 && j != 0) { 
                Point offset = Point(i * diffepsilon, j * diffepsilon);
                surroundingPoints.push_back(current + offset);
            }
            }
        }

        // cout << "current point: " << current.x << " " << current.y << endl;
        // get the gradient for each drone
        vector<double> gradient = concentration_gradient_LSQ(surroundingPoints, current);
        // cout << "gradient: " << gradient[0] << " " << gradient[1] << endl;

        double centerX = (A.x + B.x) * 0.5;
        double centerY = (A.y + B.y) * 0.5;
        double level = getGaussian(current);
        cout << "level: " << level << endl;
        cout << "contour res:" <<  contourRes << endl;
        vector<Point> localContour = getGaussianContours(level, 
            contourRes, 
            centerX - 1*DIST*epsilon, 
            centerX + 1*DIST*epsilon, 
            centerY - 1*DIST*epsilon, 
            centerY + 1*DIST*epsilon); 

        // cout << "center point: " << centerX << " " << centerY << endl;
        // cout << "DIST" << DIST << "epsilon" << epsilon << endl;
        cout << "size of local contour: " << localContour.size() << endl;

        Point criticalPoint = getCriticalPoint(localContour);
        cout << "critical point: " << criticalPoint.x << " " << criticalPoint.y << endl;

            vector<double> r = {criticalPoint.x - current.x, criticalPoint.y - current.y};
            vector<double> t = {-gradient[1], gradient[0]};

            double r_norm = sqrt(r[0] * r[0] + r[1] * r[1]);
            if (r_norm > 1e-9) {
                r[0] /= r_norm;
                r[1] /= r_norm;
            }

            double t_norm = sqrt(t[0] * t[0] + t[1] * t[1]);
            if (t_norm > 1e-9) {
                t[0] /= t_norm;
                t[1] /= t_norm;
            }


        double dotProduct = r[0] * t[0] + r[1] * t[1];
        if(dotProduct < 0){
            t[0] = -t[0];
            t[1] = -t[1];
            return t;
        }else{
            return t;
        }
        
        
    }

    std::pair<int,int> checkCross (Drone& droneA, Drone& droneB, Point motion)
    {
        // Need to get the gradient around the initial point
        // and the end point in order to compare sign of 
        // the gradient
        cout << "check cross function" << endl;
        cout << "droneA position: " << droneA.position.x + motion.x << " " << droneA.position.y + motion.y << endl;
        cout << "droneB position: " << droneB.position.x + motion.x << " " << droneB.position.y + motion.y << endl;
  

        //get tangent vector, which is normal to the gradient vector
        vector<double> tangentA = getContourTangent(droneA.position + motion, droneB.position + motion);
        vector<double> tangentB = getContourTangent(droneB.position + motion, droneA.position + motion);

        droneA.currentTangent = tangentA;
        droneB.currentTangent  = tangentB;

        droneA.currentContourGradient = {-tangentA[1], tangentA[0]};
        droneB.currentContourGradient = {tangentB[1], -tangentB[0]};

        for(int i = 0; i < tangentA.size(); i++){
            cout << "tangent A: [" << i << "]: " << tangentA[i] << endl;
        }
        
        for(int i = 0; i < tangentB.size(); i++){
            cout << "tangent B: [" << i << "]: " << tangentB[i] << endl;
        }
        

        // get direction of movement ?
        // Point droneADirection = droneA.lastMotion.getEnd() - droneA.lastMotion.getStart();
        // Point droneBDirection = droneB.lastMotion.getEnd() - droneB.lastMotion.getStart();

        // Take the dot product of the drone's motion vector with the tangent
        // Drone A is positive in counterclockwise direction (sandwhich invariant), negative in clockwise

        vector<Point> cpPoints = getCriticalPathPoints();
        vector<double> vectorBetween = {0.0, 0.0};
        if (!cpPoints.empty()) {
            Point lastCPPoint = cpPoints.back();

            // FIX ME: supposed to be vector normal to critical point,
            // but not sure how robust this is.
            vectorBetween = {-lastCPPoint.y, lastCPPoint.x}; 
        }


        double crossProductA = vectorBetween[0] * tangentA[1] - vectorBetween[1] * tangentA[0];
        double crossProductB = vectorBetween[0] * tangentB[1] - vectorBetween[1] * tangentB[0];


        cout << "cross prod A: " << crossProductA << endl;
        cout << "cross prod B: " << crossProductB << endl;

        
        std::pair<int, int> side;
        if (crossProductA > contourRes){
            side.first = RIGHT;
        }else{
            side.first = LEFT;
        }

        if (crossProductB > contourRes){
            side.second = RIGHT;
        }else{
            side.second = LEFT;
        }

        return side;


    }

    Point getCriticalPoint(vector<Point>& localContour)
    {
    
        lastCPCurvatures.clear();

        for(auto& point : localContour){
            // cout << "point: " << point.x << " " << point.y << endl;
            vector<Point> surroundingPoints;
            for(int i = -1; i < 2; i++){
                for(int j = -1; j < 2; j++){
                    if (i == 0 && j == 0) continue; 
                    Point offset = Point(i*diffepsilon, j*diffepsilon);
                    surroundingPoints.push_back(point + offset);
                }
            }

            double curvature = calc_curvature_LSQ(surroundingPoints, point); 
            // cout << "curvature: " << curvature << endl;

            lastCPCurvatures.push_back(curvature);
        }


        // get derivaties of curvature by fitting spline
        CubicSpline curv_spline;
        
        // need to get x values vs curvature
        vector<double> xValues;
        for(auto& point : localContour){
            xValues.push_back(point.x);
        }
        curv_spline.setPoints(xValues, lastCPCurvatures);
        vector<double> curv_splineDerivatives;
        for (size_t i = 0; i < xValues.size(); ++i) {
            double derivative = curv_spline.getDerivative(xValues[i]);
            curv_splineDerivatives.push_back(std::pow(derivative, 2));
            // cout << "x: " << xValues[i] << " der: " << derivative << endl;
        }

        // now get where minimized since we squared derivative
        size_t bestIdx = 0;
        double minSqVal = std::numeric_limits<double>::max();
        for (size_t i = 1; i + 1 < curv_splineDerivatives.size(); ++i) {
            // cout << "curvature derivative ** 2: " << curv_splineDerivatives[i] << endl;

            if (curv_splineDerivatives[i] < minSqVal) {
                    bestIdx = i;
                    minSqVal = curv_splineDerivatives[i];
                }
        } 



        // vector<double> curvatureDerivatives;
        // curvatureDerivatives = compute_curvature_derivative(localContour, lastCPCurvatures);

        // size_t bestIdx = 0;
        // double minSqVal = std::numeric_limits<double>::max();
        // for (size_t i = 1; i + 1 < curvatureDerivatives.size(); ++i) {
        //     // cout << "curvature derivative: " << curvatureDerivatives[i] << endl;
        //     if (std::pow(curvatureDerivatives[i],2) < minSqVal) {
        //             bestIdx = i;
        //             minSqVal = std::pow(curvatureDerivatives[i],2);
        //         }
        // }
        

        Point criticalPoint = localContour[bestIdx];

        criticalPathPoints.push_back(criticalPoint);

        return criticalPoint;
    }

    Point getCrossingPoint(Drone droneA, Drone droneB){
        // NOte that Drone A is the drone we are determining
        // the crossing point for 

        cout << "get crossing point function " << endl;
        double startLevel = getGaussian(droneA.position);
        double endLevel = getGaussian(droneA.last);
        cout << "start level: " << startLevel << endl;
        cout << "end level: " << endLevel << endl;

        // used to look for point that minimizes the curvature
        // (critical point numerically)
        auto findMinIndex = [](const std::vector<double>& v) -> size_t {
            return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
        };

        double startCenterX = (droneA.position.x + droneB.position.x) * 0.5;
        double startCenterY = (droneA.position.y + droneB.position.y) * 0.5;
        // get contour lines within box around start and end
        if (std::abs(startLevel - endLevel) > 1e-9){
            // we consider this to be two different contours
            // and will need to evaluate the critical path 
            // between two crtitical points

            // FIX ME: Need to determine size of region that we will get
            // contours for.

            double endCenterX = (droneA.last.x + droneB.last.x) * 0.5;
            double endCenterY = (droneA.last.y + droneB.last.y) * 0.5;
            
            vector<Point> localContourStart = getGaussianContours(startLevel, 
                                                             contourRes, 
                                                             startCenterX - 1*DIST*epsilon, 
                                                             startCenterX + 1*DIST*epsilon, 
                                                             startCenterY - 1*DIST*epsilon, 
                                                             startCenterY + 1*DIST*epsilon); 

            vector<Point> localContourEnd = getGaussianContours(endLevel,
                                                                contourRes,
                                                                endCenterX - 1*DIST*epsilon,
                                                                endCenterX + 1*DIST*epsilon,
                                                                endCenterY - 1*DIST*epsilon,
                                                                endCenterY + 1*DIST*epsilon);
         
            // critical path is path between two critical points                                                      
            Point criticalPointStart = getCriticalPoint(localContourStart);
            Point criticalPointEnd = getCriticalPoint(localContourEnd);

            criticalPathPoints.push_back(criticalPointStart);
            criticalPathPoints.push_back(criticalPointEnd);

            // find intersection between critical path and start and end
            double startX1 = droneA.position.getX();
            double startY1 = droneA.position.getY();
            double EndX2 = droneA.last.getX();
            double EndY2 = droneA.last.getY();
            double cpX3 = criticalPointStart.getX();
            double cpY3 = criticalPointStart.getY();
            double cpX4 = criticalPointEnd.getX();
            double cpY4 = criticalPointEnd.getY();

            double denom = (startX1 - EndX2)*(cpY3 - cpY4) - (startY1 - EndY2)*(cpX3 - cpX4);
            if (denom == 0){
                cout << "Exception! No intersection!"<<endl;
                exit(0);
            }

            double px = ((startX1 * EndY2 - startY1 * EndX2) * (cpX3 - cpX4) - (startX1 - EndX2) * (cpX3 * cpY4 - cpY3 * cpX4)) / denom;
            double py = ((startX1 * EndY2 - startY1 * EndX2) * (cpY3 - cpY4) - (startY1 - EndY2) * (cpX3 * cpY4 - cpY3 * cpX4)) / denom;
            Point intersectionPoint(px, py);

            return intersectionPoint;
                

        }else{
            // we consider this to be a single contour
            // and will need to evaluate where the critical point is
            vector<Point> localContour = getGaussianContours(startLevel, 
                                                             0.01, // FIX ME
                                                             startCenterX - 1*DIST*epsilon, //FIX ME
                                                             startCenterX + 1*DIST*epsilon, // FIX ME
                                                             startCenterY - 1*DIST*epsilon,  // FIX ME
                                                             startCenterY + 1*DIST*epsilon); //FIX ME

            Point criticalPoint = getCriticalPoint(localContour);
            criticalPathPoints.push_back(criticalPoint);
            return criticalPoint;

        }
        
    }
    
    CrossData getCross (Drone& droneA, Drone& droneB, double alpha, double dist)
    {
       /* FIX ME
        - removed inside, nabla, arg and will take it from dones passed
        - alpha and dist can stay
        - line segment data stored in each drone "lastMotion"
       */


       // FIX ME:  how to check cross for
       // each drone?
       Point motion = PointUtil::vector(droneA.nabla + alpha, dist/100);
       cout << "motion: " << motion.x << " " << motion.y << endl;
       
       for(int i = 0; i < 100; i++){

            std::pair<int,int> crossInfo = checkCross(droneA, droneB, motion);
    
            cout << "cross info" << crossInfo.first << ", " << crossInfo.second << endl;  

            Point crossPoint = Point(0,0);
            if(crossInfo.first != LEFT) // location of drone A
            {
                
                crossPoint = getCrossingPoint(droneA, droneB);
                droneA.numCross++;
                droneA.side = RIGHT;

                return CrossData (crossPoint, 1);

            }else if(crossInfo.second != RIGHT){
                
                crossPoint = getCrossingPoint(droneB, droneA);
                droneB.numCross++;
                droneB.side = LEFT;
                   
                return CrossData (crossPoint, 2);
            }
            motion += motion;
       }
       
       return CrossData ( Point (0,0), 0 );
   }

   bool foundSource(Drone droneA, Drone droneB){
        // FIX ME: if drone has encountered a concave contour we 
        // know its a source

        // FIX ME: need to check if still using this
        if (droneA.lastTangent.size() == 0){
            return false;
        }

        double crossProdLast = droneA.lastTangent[0] * droneB.lastTangent[1] - 
                               droneA.lastTangent[1] * droneB.lastTangent[0];

        double crossProdCurrent = droneA.currentTangent[0] * droneB.currentTangent[1] - 
                               droneA.currentTangent[1] * droneB.currentTangent[0];
        
        // check if the cross product of the last and current gradient
        // is negative, which indicates that the contours have
        // changed from concave to convex or vice versa
        if (crossProdLast * crossProdCurrent < 0){
            // double check if the cross product is negative
            // and positive where expected?
            cout << "found source" << endl;
            return true;
        }
        return false;


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

LineSegment::LineSegment(Point start, Point end) : line(Line::buildByPoints(start, end)), start(new Point(start.getX(), start.getY())), end(new Point(end.getX(), end.getY())) {}

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

// LineSegment::~LineSegment() {
//     delete start;
//     delete end;
// }

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

void Sync (Drone &A, Drone &B, double alpha, double dist, criticalPath &cp)
{
    cout << "syncing..." << endl;
    Point v;
    
    if (alpha > 0)
        v = PointUtil::vector (A.nabla + alpha - PI/2, DIST * epsilon);
    else
        v = PointUtil::vector (A.nabla + alpha + PI/2, DIST * epsilon);

    
    Point nextPosition = A.position + v;
    B.last = nextPosition;
    B.position = nextPosition;


    Point BLast = B.last;
    Point BPosition = B.position;

    LineSegment dronemotion = LineSegment(BLast, BPosition);   
 //   if (abs(B.nabla - A.nabla) < 1e-9 &&  plume.getCross(dronemotion, B.nabla, alpha, dist, B.inside).second == 1)
   //     B.inside = B.inside ^ 1;
    B.nabla = A.nabla;
   // B.inside = A.inside;
    B.polytope.push_back (nextPosition);
    //ignoring angle turned during sync.
    cout << "exiting sync" << endl;
    return ;
}

string check (Point P){
    double A = majorAxis;
    double B = minorAxis;
    return (((P.x * P.x) / (A * A) + (P.y * P.y) / (B*B))  <= 1) ? "Inside " : "Outside ";
}


void print_data(Drone A, Drone B)
{
    cout << "printing data..." << endl;
    int numPoints = A.polytope.size();
    
    cout << "printing data 1..." << endl;
    fprintf (out, "Pen b\n");
    
    for (int i = 0;i < numPoints; i ++)
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", A.polytope[i].x, A.polytope[i].y, A.polytope[(i+1)%numPoints].x, A.polytope[(i+1)%numPoints].y);
  
    numPoints = B.polytope.size();
  
    cout <<numPoints << endl;

//    fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[0].center.x, plume.ovals[0].center.y, R, R/2);
 //   fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[1].center.x, plume.ovals[1].center.y, R, R/2);

    cout << "printing data 2..." << endl;
    fprintf (out, "Pen r\n");
    
    for (int i = 0;i < numPoints; i ++)
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", B.polytope[i].x, B.polytope[i].y, B.polytope[(i+1)%numPoints].x, B.polytope[(i+1)%numPoints].y);
 
    cout << "done printing data..." << endl;
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

bool CrossCriticalPath(Drone &A, Drone &B, double alpha, criticalPath &cp)
{
    Point start_pos = A.last; 

    int crossing = A.side; 
    cout << "crossing " << crossing << endl;
    bool orient = true ;
    double alphainitial = alpha;
    bool endHere = false;
    
    int iterate = 1000;
    
    cout << "cross CP" << endl;

    do{
        cout << "move pair" << endl;
        cout << "alpha " << alpha << endl;
        A.MoveDrone (alpha, epsilon, cp.diffepsilon, 0);
        CrossData crossData = cp.getCross(A, B, alpha, epsilon);
        cout << "drone side? " << A.side << endl;
        cout << "found source?" << endl;
        endHere = endHere || cp.foundSource(A, B);
        
        // FIX ME: determine orientation?
        Point last = A.last;
        Point position = A.position;
        if (PointUtil::orientation(last,position, start_pos) == PointUtil::CLOCKWISE && alpha > 0)
            orient = false ;
        if (PointUtil::orientation(last, position, start_pos) == PointUtil::COUNTERCLOCKWISE && alpha < 0)
            orient = false ;

        if (alpha > 0)
            alpha += epsilon;
        else
            alpha -= epsilon;

        A.angleTurned += abs (alphainitial);

        iterate--;
        if (iterate < 0){
            print_data (A,B);
            exit (0);
        }
    }while (crossing == A.side && orient && !endHere);
        
    if (crossing == A.side && !endHere){
        A.polytope.pop_back (); // FIX ME: ?
        A.polytope.pop_back();
        A.position = A.last;
        A.angleTurned -= abs (alphainitial);
        double dx = start_pos.getX () - A.position.getX();
        double dy = start_pos.getY () - A.position.getY();
        double gradient = atan2 (dy, dx);
        
        Point pos = A.position;
        Point d1 = start_pos - pos;
        Point motion;
        if (d1.length() > epsilon*epsilon)
            motion = PointUtil::vector (gradient, epsilon * epsilon);
        else
            motion = PointUtil::vector (gradient, d1.length());
        motion = pos + motion;
        Point d2 = start_pos - motion;

        if (d1.length() < d2.length()){
            reverse (gradient);
        }

        A.angleTurned += changeGradient (A.nabla + alpha, gradient);
        A.nabla = gradient;

        int iter = 0;
        while (crossing == A.side && !endHere){
            Point new_pos = A.position; 
            d1 = start_pos - new_pos; // FIX ME: how do we want to deal with pair position
            if (d1.length() > epsilon * epsilon)
            {
                A.MoveDrone(0, epsilon*epsilon, cp.diffepsilon, 0);
                CrossData crossData = cp.getCross(A, B, alpha, epsilon);
                cout << "Drone side now " << A.side << endl;
                endHere = endHere || cp.foundSource(A, B); 
            }
            else{
                A.MoveDrone(0, d1.length(), cp.diffepsilon, 0);
                CrossData crossData = cp.getCross(A, B, alpha, epsilon);
                cout << "Drone side now " << A.side << endl;
                endHere = endHere || cp.foundSource(A, B);
            }
            
            ++iter;
            if (iter > 10000)
            {
                cout<<"Iterations exceeding ..."<<endl;
                print_data(A,B);
                exit (0);
            }
        }
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
    cout << "Running Sketch Algorithm for eapsilon = " << epsilon << endl;
    
    vector<Ellipse> ell;
    
    for (int i = 0;i < num; i ++)
        ell.push_back (Ellipse(gaussianCenter[i], majorAxis, minorAxis));
    
    PLUME plume; //FIX ME: delete
    
    
    /*
     FIX ME
     
     INITIALIZE THE ALGORITHM WITH DRONE PAIRS
     AND APPROPRIATE LOCATION OF THE PAIRS
    
     */
    
    plume.ovals = ell; 
   // fprintf (out, "Ellipse (%lf,%lf) %lf %lf ", plume.ovals[0].center.x, plume.ovals[0].center.y, majorAxis, minorAxis);
    
    Drone A (drone_start_A, drone_start_A, 1, 1.47, true);
    // Drone ab (drone_start_AB, drone_start_AB, 1, PI/2, true);

    Drone B (drone_start_B, drone_start_B, 2, 1.47, false);
    // Drone bb (drone_start_BB, drone_start_BB, 2, PI/2, false);

    //create drone pairs 
    // FIX ME: set lambda
    // DronePair A(aa, ab, 0); 
    // DronePair B(ba, bb, 0);
    
    int TOKEN = 2;
    
    bool loopEnd = false ;

    // FIX ME: remove start point - uneeded
    Point startPoint = {(drone_start_A.x + drone_start_B.x) / 2, (drone_start_A.y + drone_start_B.y) / 2};
    double res = epsilon*epsilon; // how much resolution should we calculate contourlines at?
    criticalPath cp(startPoint, epsilon, res);

    alpha = 0; // test
    // CrossData crossData = cp.getCross(A, B, alpha, epsilon); // just to initialize tangent and graidnet data

    
  //  fprintf (out, "Pen b\n");
    
    
    do{
        cout << "top of do while " << endl;
        int iter = 0;
        cout << "iteration " << iter << endl;

       // Ensure A and B are properly initialized before entering the loop

        //FIX ME: CHANGE A, B TO DRONE PAIRS AND IN THE FUNCTION CALL CROSS PLUME.
        while ((A.numCross + B.numCross == 0 || 3 == A.side + B.side) && !loopEnd)
        {
         //   cout << "testing TOKEN ... "<< A.position.x << " "<< A.position.y <<" 1"<<  endl;
            cout << "inside loop " << iter << endl;
            ++iter;
            if (iter > 10000)
            {
                cout<<"Iterations exceeding ..."<<endl;
                print_data(A,B); // FIX ME
                exit (0);
            }
            
            // FIX ME: get terminating condition from collection of
            // values passed from Move Drone
            // loopEnd = loopEnd || A.MovePair (alpha, epsilon, plume, 1);
            // loopEnd = loopEnd || B.MovePair (alpha, epsilon, plume, 1);

            cout << "moving drones " << endl;
            cout << "alpha " << alpha << endl;
            A.MoveDrone(alpha, epsilon, cp.diffepsilon, 1);
            B.MoveDrone(alpha, epsilon, cp.diffepsilon, 1);
            cout << "moved...now see if they found source" << endl;
            loopEnd = loopEnd || cp.foundSource(A, B);
            
            // FIX ME: Modify MoveDrone to return gradient contour
            // and determine inside/outside from the returned gradient
            // We will call the getCross function from the criticalPath
            // and pass the info about each drone to get whether or
            // not they drones crossed, if they are inside or outside
            cout << "checking crossing ..." << endl;
            CrossData crossData = cp.getCross(A, B, alpha, epsilon);
            cout << "cross data " << crossData.first.getX() << ", " << crossData.second << endl;
            cout << "finished checking crossing" << endl;
            cout << "crossData.first " << crossData.first.getX() << ", " << crossData.first.getY() << endl; 
            cout << "crossData.second " << crossData.second << endl;
            if(crossData.second){
                cout << "numcross A" << A.numCross << endl;
                cout << "numcross B" << B.numCross << endl;
                cout << "A side" << A.side << endl;
                cout << "B side" << B.side << endl;
                // call function to get gradient at crossing point
                // FIX ME: need to determine which drone crossed
                // and then get the gradient
                cout << "crossed; learning gradient..." << endl;
                vector<double> gradient_vec = cp.getGradientAtPoint(crossData.first);
                A.LearnGradient(alpha, epsilon, crossData.first, B, gradient_vec);
                B.LearnGradient(alpha, epsilon, crossData.first, A, gradient_vec);

            }

            cout << "endloop" << endl;
            
        }

        cout << "left inner loop " << endl;

        if (A.side + B.side != 3)
        {
            cout << "A.getInside() + B.getInside() != 3" << endl;

            // If A crosses
            if (A.side == 2 && B.side == 2)
            {
                cout << "A.getInside() & B.getInside() == 2" << endl; 
                alpha = -epsilon;
            //    B.nabla = A.nabla;

                loopEnd = loopEnd || CrossCriticalPath (B,A, alpha, cp);
              //  if (A.inside + B.inside == 1)
                cout << "about to sync" << endl;
                Sync (B,A,alpha, epsilon, cp); // FIX ME
            
                    A.nabla = B.nabla;
              
 //               cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
 //               cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
            }
            else // B crosses
            {
                cout << "inner if else" << endl;
            //    cout<<"Exception both inside!"<<endl;
          //      exit(0);
                alpha = epsilon;
        //        Sync (B,A,alpha, epsilon, plume);
              //  A.nabla = B.nabla;

                
                loopEnd = loopEnd || CrossCriticalPath (A,B, alpha, cp);
              //  if (A.inside + B.inside == 1)
              cout << "about to sync" << endl;
                Sync (A,B,alpha, epsilon, cp);

                cout << "set B nabla " << endl;
                    B.nabla = A.nabla;
                
   //             cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
   //             cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
            
            }
            cout << "exit if" << endl;
        }
        cout << "exit outer if" << endl;
    }while (!loopEnd);
    
    cout << "exited loop" << endl;

    // FIX ME
    if (!A.polytope.empty()) {
        cout <<"Initial crossing with A is " << A.polytope[0].getX() << " " << A.polytope[0].getY() << endl;
    } else {
        cout << "Initial crossing with A is not available as polytope is empty." << endl;
    }
    cout << "angle turned by A is " << " " << A.angleTurned << endl;
    cout << "distance traversed by A is "<< " " << A.distTraversed << endl;
    if (!A.polytope.empty()) {
        cout << "area estimated by A is " << " " << estimateArea (A.polytope) << endl;
        areas.push_back (estimateArea (A.polytope));
    } else {
        cout << "area estimated by A is not available as polytope is empty." << endl;

    lengths.push_back (A.distTraversed);
    angles.push_back (A.angleTurned);
    cout << "actual area is  " << PI * majorAxis * minorAxis << endl;
    }

    print_data (A,B);
    
    return ;
}   


void test_infrastructure()
{
    gaussianCenter.clear();
    gaussianVar.clear();
    double alpha = 0; // PI/2;
    num = 4;
    CROSSBOUND = 100;
    majorAxis = 0.25;
    minorAxis = 0.25;
    DIST = sqrt (49);
    THRESHOLD = exp (-majorAxis*majorAxis);
   
    // 0.9225 -3.39337
    drone_start_B = Point (1 + DIST*epsilon*0.5,-2.38);
    drone_start_A = Point (1 - DIST*epsilon*0.5,-2.38);

    drone_start_BB = Point (1 + DIST*epsilon*0.5,-2.38);
    drone_start_AB = Point (1 - DIST*epsilon*0.5,-2.38);
    
    gaussianCenter.push_back(Point(1,1));
    gaussianCenter.push_back(Point(-2,0));
    gaussianCenter.push_back(Point(1,0));
    gaussianCenter.push_back(Point(4,0));
    
    gaussianVar.push_back(Point(0.6,0.6));
    gaussianVar.push_back(Point(0.8,1.0));
    gaussianVar.push_back(Point(0.9,1.8));
    gaussianVar.push_back(Point(0.5,0.5));

    cout<<"Initializing test... "<<endl;
    cout<<"Epsilon "<<epsilon<<endl;
    cout<<"Initial direction "<<alpha<<endl;
    cout<<"Number of Gaussians "<<num<<endl;
    cout<<"Least difference between starting and end point "<<INF<<endl;
    cout<<"Minimum crossings before checking termination "<<CROSSBOUND<<endl;
   // cout<<"Major and Minor Axis of Gaussian "<<majorAxis<<" "<<minorAxis<<endl;
    cout<<"Minimum distance factor between drones and minimum distance "<<DIST<<" " <<DIST * epsilon<<endl;
    cout<<"Concentration THRESHOLD "<<THRESHOLD<<endl;
    cout<<"Drone A starting point "<<drone_start_A.x<<" "<<drone_start_B.y<<endl;
    cout<<"Drone B starting point "<<drone_start_B.x<<" "<<drone_start_A.y<<endl;
    
    cout<<"Centers of gaussians "<<endl;
    
    for (int i =0 ; i < num; i ++)
        cout<<gaussianCenter[i].x<<" "<<gaussianCenter[i].y<<endl;
        
    sketch_algorithm(alpha);
}

int main()
{
    int i = 0;
    // 0.005 originally
    for (epsilon = 0.01; i < 1 ; epsilon += 0.001){
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
