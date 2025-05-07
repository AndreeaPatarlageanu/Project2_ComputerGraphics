//what I write from here to the class Polygon will be written in class
//to be rearranged after the class
//to add the Vector .h and .cpp from the previous project!!!

#include <random>
#include <vector>
#include "Vector.h"
#define M_PI 3.14159265359

static std::default_random_engine engine(10);  //random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

const double EPSILON = 1e-4;

void boxMuller( double stdev, double &x, double &y ){
	double r1 = uniform( engine );
	double r2 = uniform( engine );
	x = sqrt( -2 * log(r1)) * cos( 2 * M_PI * r2 ) * stdev;
	y = sqrt( -2 * log(r1)) * sin( 2 * M_PI * r2 ) * stdev;
}

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector operator*( const Vector& a, const Vector& b ) {
    return Vector( a[0] * b[0], a[1] * b[1], a[2] * b[2] );
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    std::vector<Vector> vertices;
};	


//add all the vector studd here
class VoronoiDiagram{
	public:
	VoronoiDiagram() {};

	Polygon clip_by_bisector( const Polygon& V, const Vector&P0, const Vector& Pi ){
		Polygon result;
		const Vector P0Pi = Pi - P0;
		for( int i = 0; i < V.vertices.size(); i++ ){

			const Vector &A = V.vertices[ ( i == 0 ) ? V.vertices.size() -1: i - 1 ];
			const Vector &B = V.vertices[ i ];

			if( ( B - P0 ).norm2() <= ( B - Pi ).norm2() ){  //B is inside
				if( ( A - P0 ).norm2() >= ( A - Pi ).norm2() ){
					Vector M = ( P0 + Pi ) * 0.5;
					double t = dot( M - A, P0Pi) / dot( B - A, P0Pi );
					Vector P = A + t * ( B - A );
					result.vertices.push_back( P );
				}
				result.vertices.push_back( B );
			}
			else {
				if( ( A - P0 ).norm2() <= ( A - Pi ).norm2() ){
					Vector M = ( P0 + Pi ) * 0.5;
					double t = dot( M - A, P0Pi ) / dot( B - A, P0Pi );
					Vector P = A + t * ( B - A );
					result.vertices.push_back( P );
				}
			}
		}
		return result;
	}

	void compute() {
		Polygon square;

		square.vertices.push_back(Vector(0, 0));
		square.vertices.push_back(Vector(0, 1));
		square.vertices.push_back(Vector(1, 1));
		square.vertices.push_back(Vector(1, 0));

		diagram.resize(points.size());

        #pragma omp parallel for schedule( dynamic, 1)
		for( int i = 0; i < points.size(); i++ ) {
			Polygon V = square;
			for ( int j = 0; j < points.size(); j++ ){
				if( i == j ) continue;

				V = clip_by_bisector(V, points[i], points[j]);
			}
			diagram[i] = V;
		}
	}
	std::vector<Vector> points;
	std::vector<Polygon> diagram;

};


// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points = NULL, std::string fillcol = "none") {
	FILE* f = fopen(filename.c_str(), "w+");
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (int i = 0; i < polygons.size(); i++) {
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \"");
		for (int j = 0; j < polygons[i].vertices.size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
	}

	if (points) {
		fprintf(f, "<g>\n");		
		for (int i = 0; i < points->size(); i++) {
			fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"3\" />\n", (*points)[i][0]*1000., 1000.-(*points)[i][1]*1000);
		}
		fprintf(f, "</g>\n");

	}

	fprintf(f, "</svg>\n");
	fclose(f);
}

int main() {
	int N = 100;

	VoronoiDiagram Vor;
	
	for( int i = 0; i < N; i++ ) {
		Vor.points.push_back( Vector(uniform(engine), uniform(engine), 0.0 ) );
	}

	Vor.compute();

	save_svg( Vor.diagram, "testOut2.svg", &Vor.points);

}
