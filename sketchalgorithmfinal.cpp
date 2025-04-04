#include <cmath>
#include <vector>
#include <cstdio>
#include <random>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

const int ITER = 10;
const int CROSSBOUND = 15;
double MAXVAL = (1e10) + 1;
const double R = 0.5;
const double PI = 3.14159;
const double DIST = sqrt (2);
const double X [6] = {0,1,0.75,1,0,0.25};
const double Y [6] = {0,0,0.5,1,1,0.5};

double epsilon, INF;

vector<double> areas, lengths, angles, eps;

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

typedef pair<Point, int> CrossData;

double get_dist (Point A, Point B)
{
    return sqrt ((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
}



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

struct LineSegment {

    LineSegment(Point &start, Point &end);

    LineSegment(const Line &line, const Point &start, const Point &end);

    LineSegment(const LineSegment& copySegment);

    ~LineSegment();
    
    bool inside (Point P)
    {
        double d = length();
        double d1 = get_dist (*start, P);
        double d2 = get_dist (P, *end);
        if (abs (d1 + d2 - d) < 1e-9)
            return true ;
        return false;
    }
    
    Point intersection (LineSegment S) //the argument here is always the dronePath
    {
        Line L = getLine();
        Line M = S.getLine();
        double x, y;
        
        if (L.getM() > 1e10 && M.getM() < 1e10){
            x = L.getB();
            y = x * M.getM() + M.getB();
            return Point (x,y);
        }
        
        if (abs(L.getM() - M.getM()) < 1e-9)
        {
            if (abs(L.getB() - M.getB()) < 1e-9){
                if (S.inside (*start))
                    return *start ;
                if (S.inside (*end))
                    return *end;
                if (inside (*S.start))
                    return *S.start;
                if (inside(*S.end))
                    return *S.end;
                else
                    return Point (100,100);
            }
            else
                return Point (100,100);
        }
        
        x = (L.getB()-M.getB())/(M.getM()-L.getM());
        y = L.getM() * x + L.getB() ;
        
        return Point (x,y);
    }
    
    double getGradient ()
    {
        double gradient = atan2(end->y - start->y, end->x - start->x);
        Point motion = PointUtil::vector (gradient, length());
        Point next = *start + motion;
        
        if (abs(next.y - end->y) > 1e-9 || abs (next.x - end->x) > 1e-9)
            reverse (gradient);
        
        return gradient;
    }
    
    double length();

    Line getLine();

    Point getStart();

    Point* getStartPtr();

    Point getEnd();

    Point* getEndPtr();

    Line line;
    Point *start, *end;
};

double changeGradient (double angleA, double angleB)
{
    return abs (angleA - angleB); //check this again for bug.
}

struct Graph
{
    vector<Point> nodes;
    vector<LineSegment> edges;
    Graph()
    {
        
    }
    Graph(vector<Point> V)
    {
        nodes = V;
        int n = V.size();
        for (int i = 0;i < V.size() ; i ++)
            edges.push_back (LineSegment (V[i], V[(i+1)%n]));
    }
    bool intersect (LineSegment A, LineSegment B)
    {
        Point common = A.intersection (B);
   
       // if (abs(common.x) < 1 + epsilon && abs(common.y) < 1 + epsilon)
       // cout<<"Common " << common.x<<" "<<common.y<<" "<<B.inside(common)<<" "<<" "<<B.start->x<<" "<<B.start->y<<" "<<B.end->x<<" "<<B.end->y<<" "<<B.length()<<endl;
        
        if (A.inside (common) && B.inside(common))
            return true ;
        return false ;
    }
    bool crosses (LineSegment DronePath)
    {
        for (int i = 0;i < edges.size(); i++)
            if (intersect (edges[i], DronePath))
                return true ;
        return false;
    }
    LineSegment getCross(LineSegment DronePath)
    {
        for (int i = 0;i < edges.size(); i ++)
        {
            if (intersect(edges[i], DronePath)){
         //       cout<<"Edge is "<<i<< " "<<DronePath.start->x<<" "<<DronePath.start->y<<" "<<DronePath.end->x<<" "<<DronePath.end->y<<endl;
           //     cout<<"Edge data "<<edges[i].start->x<<" "<<edges[i].start->y<<endl;
                return edges[i] ;
            }
        }
        cout << "Does not intersect, exception! " << endl;
        return edges[0];
    }
};


class Drone {
    public :
    Point position;
    Point last;
    double nabla ; //just the gradient angle, not the slope
    int inside, numCross ;
    bool droneIn;
    double angleTurned;
    double distTraversed;
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
    }
    
    bool MoveDrone (double alpha, double dist, Graph &plume)
    {
      //  cout<<"moving drone .. "<<endl;
        Point motion;
        Point nextPosition;
        motion = PointUtil::vector (nabla + alpha, dist);
        nextPosition = position + motion;
        LineSegment dronemotion = LineSegment (position, nextPosition);
        
        distTraversed += dist;
        
        if (abs(nextPosition.x) > 3 || abs(nextPosition.y) > 3)
            cout << "position exception! "<<endl;
        
        Point crossingPoint;
        int cross = plume.crosses (dronemotion);
        if (cross)
        {
            ++numCross;
            inside = inside ^ 1;
            LineSegment crossingsegment = plume.getCross (dronemotion);
            double gradient = crossingsegment.getGradient();
            crossingPoint = crossingsegment.intersection(dronemotion);
            if (droneIn)
            {
                cout<<"Crossed at "<<crossingPoint.x<< " "<< crossingPoint.y<<" " <<gradient<<" "<<droneIn<<endl;
                cout<<"NEXT POSITION "<<nextPosition.x<<" "<<nextPosition.y<<endl;
            }
            angleTurned += changeGradient (nabla + alpha, gradient);
            nabla = gradient;
        }
        last = position;
        position = nextPosition;
        if (polytope.size())
            polytope.push_back (position);
        if (numCross == 1 && cross)
            polytope.push_back (crossingPoint);
        
      //  if (get_dist (position, last) < 1e-9)
     //       cout << "Exception with distance!"<<endl;
        
        if (polytope.size() > 1)
        {
            Point currtoinit = polytope.back() - polytope[0];
         //   if (numCross >= CROSSBOUND && droneIn)
           // cout<<currtoinit.length() << " loop end "<<polytope[0].x<<" "<<polytope[0].y<<" "<<polytope.back().x<<" "<<polytope.back().y<<endl;
           // cout << currtoinit.length()  << " " << numCross << endl;
            return (currtoinit.length() < INF) && (numCross >= CROSSBOUND);
        }
        else
            return false ;
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
    if (abs(end.getX() - start.getX()) < 1e-9){
        return Line (MAXVAL, start.getX());
    }
    double m = (end.getY() - start.getY()) / (end.getX() - start.getX()); //divide by zero case solved by 1e-9
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

void Sync (Drone &A, Drone &B, double alpha)
{
    Point v;
    
    if (alpha > 0)
        v = PointUtil::vector (A.nabla + alpha - PI/2, DIST * epsilon);
    else
        v = PointUtil::vector (A.nabla + alpha + PI/2, DIST * epsilon);
    
    Point nextPosition = A.position + v;
    B.last = nextPosition;
    B.position = nextPosition;
    B.nabla = A.nabla;
    B.polytope.push_back (nextPosition);
    //ignoring angle turned during sync.
    
    return ;
}
/*
string check (Point P){
    double A = R;
    double B = R/2;
    return (((P.x * P.x) / (A * A) + (P.y * P.y) / (B*B))  <= 1) ? "Inside " : "Outside ";
}*/

bool CrossPlume (Drone &A, Drone &B, double alpha, Graph &plume)
{
    Point start_pos = A.last;
   // if (abs(start_pos.x) < 2 && abs(start_pos.y) < 2)
     //   cout <<"start pos ... "<< start_pos.x << " " << start_pos.y <<" "<<check(A.position)<< endl;

    int crossing = A.inside;
    bool orient = true ;
    double alphainitial = alpha;
    bool endHere = false;
    
    do{
        endHere = endHere || A.MoveDrone (alpha, epsilon * epsilon, plume);
        
        if (PointUtil::orientation(A.last, A.position, start_pos) == PointUtil::CLOCKWISE && alpha > 0)
            orient = false ;
        if (PointUtil::orientation(A.last, A.position, start_pos) == PointUtil::COUNTERCLOCKWISE && alpha < 0)
            orient = false ;
        Sync (A,B,alpha);
        if (alpha > 0)
            alpha += epsilon;
        else
            alpha -= epsilon;
        A.angleTurned += abs (alphainitial);
    //    cout << "testing cross plume ... "<< A.position.x << " " << A.position.y <<" "<<orient<<" "<<A.inside<<" "<<crossing<<" "<<endHere<< endl;
       // cout << "testing cross plume ... "<< B.position.x << " " << B.position.y << endl;

    }while (crossing == A.inside && orient && !endHere);
 
  //  cout << "testing cross plume exit ... "<< A.position.x << " " << A.position.y <<" "<<orient<<" "<<A.inside<<" "<<crossing<<" "<<endHere<< endl;

    if (crossing == A.inside && !endHere){
        A.polytope.pop_back ();
        cout <<"Polytope size "<< A.polytope.size() << endl;
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
       // cout<<A.position.x<<" " <<A.position.y<<" "<<d1.length()<< " distances " << d2.length()<<endl;

        if (d1.length() < d2.length()){
            reverse (gradient);
        }
        A.angleTurned += changeGradient (A.nabla + alpha, gradient);
        A.nabla = gradient;
     //   B.nabla = gradient;
     //   cout << start_pos.x << " starting here " << start_pos.y << " "<<gradient<<" "<<check (start_pos) <<" "<<check(A.position)<<" "<<A.inside<< endl;
     //   cout << "testing cross plume2 ... "<< A.position.x << " " << A.position.y << endl;
     //   cout << "testing cross plume2 ... "<< A.last.x << " " << A.last.y << endl;

        while (crossing == A.inside && !endHere){
            d1 = start_pos - A.position;
            if (d1.length() > epsilon * epsilon)
            {
                endHere = endHere || A.MoveDrone (0, epsilon*epsilon, plume);
             //   B.MoveDrone (0, epsilon*epsilon, plume);
            }
            else{
                endHere = endHere || A.MoveDrone (0, d1.length(), plume);
            //    B.MoveDrone (0, d1.length()/2, plume);
            }
            
            if (abs(A.position.x) > 2 || abs(A.position.y) > 2)
                break ;
            
        //    Sync2 (A,B, alpha);
       //     cout << "testing cross plume2 ... "<< A.position.x << " " << A.position.y << endl;
       //     cout << "testing cross plume2 ... "<< B.position.x << " " << B.position.y << endl;
        }
      //  B.nabla = A.nabla;
        Sync (A,B,alphainitial);
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

void normalize (vector<Point> &input)
{
    double minX, minY;
    double maxX, maxY;
    
    minX = 1e9;
    minY = 1e9;
    maxX = 0;
    maxY = 0;
    
    double radii;
    
    for (int i = 0; i < input.size(); i ++)
    {
        minX = min (minX, input[i].x);
        maxX = max (maxX, input[i].x);
        minY = min (minY, input[i].y);
        maxY = max (maxY, input[i].y);
    }
    
    radii = max (maxX - minX, maxY - minY);
    
    for (int i = 0;i < input.size(); i ++)
    {
        input[i].x -= minX;
        input[i].y -= minY;
        input[i].x /= radii;
        input[i].y /= radii;
    }
    
    return ;
}

void sketch_algorithm ()
{
    FILE *in = fopen ("sketchinput.txt", "r");
    FILE *out = fopen ("sketch_plot.txt", "w");
    
    cout << "Running Sketch Algorithm for Epsilon = " << epsilon << endl;
    
    vector<Point> input ;
    
  //  int SIDES = 100;
    
    // uncomment this part to activate regular polygon input
    //for (int i = 0; i < SIDES; i ++){
   //     input.push_back (Point (cos(i*(2*PI/SIDES))/2, sin(i*(2*PI/SIDES))/2));
   // }
    
    int N ;
    
    //comment this part to deactivate concave polygon
 //   for (int i = 0 ;i < N; i ++)
   //     input.push_back (Point(X[i], Y[i]));
    
    fscanf (in, "%d", &N);
    
    for (int i = 0;i < N; i ++)
    {
        int x , y ;
        fscanf (in, "%d %d", &x, &y);
        input.push_back (Point (x,y));
        cout << x << "  " << y << endl;
    }
    
    normalize (input);
    
    for (int i = 0;i < N; i++)
        cout << input[i].x<< " "<< input[i].y<<endl;
    
    Graph plume (input);
    
    double startX = 0.6;
    double startY = 0.15;
    
    //starting point for concave polygon
    Drone A (Point (startX, startY), Point (startX,startY), 1, 0, true);
    Drone B (Point (startX,DIST * epsilon + startY), Point (startX,DIST * epsilon + startY), 1, 0, false);
    
    
    //starting point for regular polygon
   // Drone B (Point (0,epsilon), Point (0,epsilon ), 1, 0, false);
  //  Drone A (Point (0,(DIST +1)* epsilon), Point ( 0, (DIST + 1)*epsilon), 1, 0, true);
     
    
 /*   fprintf (out, "Pen r\n");
    
    for (int i = 0;i < input.size(); i ++)
    {
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", input[i].x, input[i].y, input[(i+1)%N].x, input[(i+1)%N].y);
    }
   */
    int TOKEN = 2;
    
    double alpha = 0;
    
    bool loopEnd = false ;
    
    do{
        while ((A.numCross + B.numCross == 0 || 1 == A.inside + B.inside) && !loopEnd )
        {
         //   cout << "testing TOKEN A... "<< A.position.x << " "<< A.position.y <<" "<<A.inside<<" "<< B.inside<<" "<<A.nabla<<" "<<alpha<<  endl;
         //   cout << "testing TOKEN B... "<< B.position.x << " "<< B.position.y <<" 1"<<  endl;


          //  if (abs(A.position.x) > 2 || abs(A.position.y) > 2)
            //    break ;
            
         //   if (abs(B.position.x) > 2 || abs(B.position.y) > 2)
           //     break ;
          //  if (loopEnd) cout<<"Loop ended!"<<endl;
            loopEnd =  A.MoveDrone (alpha, epsilon, plume) || loopEnd;
         //   Sync (A,B,alpha);
            loopEnd =  B.MoveDrone (alpha, epsilon, plume) || loopEnd;
            
        }
        if (A.inside + B.inside != 1)
        {
            if (A.inside + B.inside == 0)
            {
              //  cout<<A.numCross <<  " number of crossings " << B.numCross<<endl;
                alpha = -epsilon;
                loopEnd = loopEnd || CrossPlume (A,B, alpha, plume);
                B.nabla = A.nabla;
                //cout<<A.inside <<  " number of insides " << B.inside<<endl;


            //    cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
            //    cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
            }
            else
            {
            //    cout<<A.numCross <<  " number of crossings " << B.numCross<<endl;

                alpha = epsilon;
                loopEnd = loopEnd || CrossPlume (B,A, alpha, plume);
                A.nabla = B.nabla;
                
              //  cout<<A.inside <<  " number of insides " << B.inside<<endl;


        //        cout << "testing Sync B... "<< B.position.x << " "<<B.position.y <<" "<<alpha<<" "<<B.nabla<< endl;
        //        cout << "testing Sync A... "<< A.position.x << " "<<A.position.y <<" "<<alpha<<" "<<A.nabla<< endl;
            }
        }
    //    TOKEN = A.inside + B.inside;
     //   lineA = A.polytope.back() - A.polytope[0];
     //   lineB = B.polytope.back() - B.polytope[0];
      //  cout << "position and nabla "<<A.position.x << " "<< A.position.y << " "<<A.nabla<< endl;
     //   if (A.polytope.back().x < 2 && A.polytope.back().y < 2 )
    //    cout << "position of polytope "<<A.polytope.back().x << " "<<A.polytope.back().y << endl;
     //   cout << lineA.length() << " distance from origin" << endl;
    }while (!loopEnd);

    cout <<"Initial crossing with A is " << A.polytope[0].x << " " << A.polytope[0].y<< endl;
    cout << "angle turned by the algorithm is " << " " << A.angleTurned + B.angleTurned  << endl;
    cout << "distance traversed by A is "<< " " << A.distTraversed << endl;
    cout << "area estimated by A is " << " " << estimateArea (A.polytope) << endl;
    areas.push_back (estimateArea (A.polytope));
    lengths.push_back (A.distTraversed);
    angles.push_back (A.angleTurned  + B.angleTurned);
    cout << "actual area is pending calculation "  << endl;
    
    int numPoints = B.polytope.size();
    
    cout <<numPoints << endl;
  
//    fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[0].center.x, plume.ovals[0].center.y, R, R/2);
 //   fprintf (out, "Ellipse (%lf,%lf) %lf %lf \n", plume.ovals[1].center.x, plume.ovals[1].center.y, R, R/2);
 
    fprintf (out, "Pen r\n");
    
    for (int i = 0;i < B.polytope.size(); i ++)
    {
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", B.polytope[i].x, B.polytope[i].y, B.polytope[(i+1)%numPoints].x, B.polytope[(i+1)%numPoints].y);
    }
    
    numPoints = A.polytope.size();
    
    fprintf (out, "Pen b\n");
    
    for (int i = 0;i < A.polytope.size(); i ++)
    {
        fprintf (out, "Line (%lf,%lf) (%lf,%lf)\n", A.polytope[i].x, A.polytope[i].y, A.polytope[(i+1)%numPoints].x, A.polytope[(i+1)%numPoints].y);
    }
   
    return ;
}

int main()
{
    int i = 0;
    for (epsilon = 0.01; i < 1 ; epsilon += 0.001){
        INF = 4 * epsilon;
     //   if (epsilon == 0.002) continue;
        sketch_algorithm ();
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
