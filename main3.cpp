//what I write from here to the class Polygon will be written in class
//to be rearranged after the class
//to add the Vector .h and .cpp from the previous project!!!
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include "Vector.h"
#include "lbfgs.h"
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
    double area() {
        if( vertices.size() < 3 )
            return 0;

        double s = 0.0;

        for( int i = 0; i < vertices.size(); i++ ) {
            int ip = ( i == vertices.size() - 1 ) ? 0 : ( i + 1 );
            //s += std::abs( vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1] );
            s += vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1];
        }
        return std::abs( s ) / 2;
    }

    double integral_square_distance(const Vector& Pi ){
        if( vertices.size() < 3 ) 
            return 0;

        double s = 0;

        for( int t = 1; t < vertices.size() - 1; t++ ){
            Vector c[3] = { vertices[0], vertices[t], vertices[t + 1] };

            double integral = 0;

            for( int k = 0; k < 3; k++ ){
                for( int l = k; l < 3; l++ ){
                    integral += dot( c[k] - Pi, c[l] - Pi );
                }
            }
            Vector edge1 = c[1] - c[0];
            Vector edge2 = c[2] - c[0];
            //double areaT = (c[1][1] - c[0][1]) *( c[2][0]-c[0][0] ) ..(POZA)
            double areaT = 0.5 * std::abs( edge1[0] * edge2[1] - edge1[1] * edge2[0] );
            s += integral * areaT / 6.;
        }
        return s;
    }
    
};	


//add all the vector studd here
class VoronoiDiagram{
public:
	VoronoiDiagram() {};

	Polygon clip_by_bisector( const Polygon& V, const Vector&P0, const Vector& Pi, const double w0, const double wi ){
		Polygon result;
		const Vector P0Pi = Pi - P0;
        double two_sqrdit = 2 * P0Pi.norm2();

        Vector M = ( P0 + Pi ) * 0.5;
        Vector Mprime = M + ( w0 - wi ) / ( two_sqrdit ) * P0Pi;
        

		for( int i = 0; i < V.vertices.size(); i++ ){

			const Vector &A = V.vertices[ ( i == 0 ) ? V.vertices.size() -1: i - 1 ];
			const Vector &B = V.vertices[ i ];

			if( ( B - P0 ).norm2() - w0 <= ( B - Pi ).norm2() - wi ){  //B is inside

				if( ( A - P0 ).norm2() - w0 >= ( A - Pi ).norm2() - wi){
					
					double t = dot( Mprime - A, P0Pi ) / dot( B - A, P0Pi );
					Vector P = A + t * ( B - A );
					result.vertices.push_back( P );
				}

				result.vertices.push_back( B );
			}
			else {

				if( ( A - P0 ).norm2() - w0 <= ( A - Pi ).norm2() - wi){
					//Vector M = ( P0 + Pi ) * 0.5;
					double t = dot( Mprime - A, P0Pi ) / dot( B - A, P0Pi );
					Vector P = A + t * ( B - A );
					result.vertices.push_back( P );
				}

			}
		}
		return result;
	}

	void compute() {
		Polygon square;

		square.vertices.push_back( Vector( 0, 0 ) );
		square.vertices.push_back( Vector( 0, 1 ) );
		square.vertices.push_back( Vector( 1, 1 ) );
		square.vertices.push_back( Vector( 1, 0 ) );

		diagram.resize( points.size() );

        #pragma omp parallel for schedule( dynamic, 1 )
		for( int i = 0; i < points.size(); i++ ) {
			Polygon V = square;
			for ( int j = 0; j < points.size(); j++ ){
				if( i == j ) continue;

				V = clip_by_bisector( V, points[i], points[j], weights[i], weights[j] );
			}
			diagram[i] = V;
		}
	}
	std::vector<Vector> points;
    std::vector<double> weights;
	std::vector<Polygon> diagram;

};

class OptimalTransport{
public:
    OptimalTransport(){};
    VoronoiDiagram vor;
    void optimize(){
        int N = vor.weights.size();
        lbfgsfloatval_t fx;
        std::vector<double> weights( N, 0 );

        lbfgs_parameter_t param;
        lbfgs_parameter_init( &param );

        int ret = lbfgs( N, &weights[0], &fx, _evaluate, _progress, ( void*)this, &param );
        memcpy( &vor.weights[0], &weights[0], N * sizeof( weights[0] ) );
        //vor.compute();
    }

    //aici bagam fct din sample.c care trb downloadat de pe git
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OptimalTransport*>( instance )->evaluate( x, g, n, step );
    }

    //DELETE STATIC!
    lbfgsfloatval_t evaluate(
        //void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        //memcpy(&ot->vor.weights[0], &weights[0], N * sizeof(weights[0]));
        //ot->vor.compute();
        memcpy( &( vor.weights[0] ), x, n * sizeof( x[0] ) );
        vor.compute();
        lbfgsfloatval_t fx = 0.0;
        for( int i = 0; i < n; i++ ) {
            double current_area = vor.diagram[i].area();
            g[i] = -( 1. / n - current_area );

            fx += vor.diagram[i].integral_square_distance(vor.points[i]) - x[i] * ( current_area - 1.0 / n);

        }
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OptimalTransport*>( instance )->progress( x, g, fx, xnorm, gnorm, step, n, k, ls );
    }
    //DELETE STATIC
    int progress(
        //void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }
};

//i guesshere i add the other classes? idk i am so lost

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
	int N = 1000;

	VoronoiDiagram Vor;
	
	for( int i = 0; i < N; i++ ) {
		Vor.points.push_back( Vector( uniform( engine ), uniform( engine ), 0.0 ) );
        Vor.weights.push_back( 0 );
	}

    auto start = std::chrono::high_resolution_clock::now();

	Vor.compute();

    OptimalTransport ot;
    ot.vor = Vor;
    ot.optimize();

    auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	std::cout << "Voronoi compute time: " << elapsed.count() << " seconds" << std::endl;

	save_svg(ot.vor.diagram, "testOut_equal_1000.svg", &ot.vor.points);

}
