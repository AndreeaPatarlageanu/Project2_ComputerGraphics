#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "Vector.h"

#define RAND_MAX() 2147483647

//I was helped by Martin Lau for this ungraded TD

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

void get_random_direction() {
	return Vector(
		((double) rand() / RAND_MAX()) * 2 - 1,
		((double) rand() / RAND_MAX()) * 2 - 1,
		((double) rand() / RAND_MAX()) * 2 - 1
	);
}

int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	//the other image as well
	unsigned char *image2 = stbi_load( "redim.jpg", 
								&W,
								&H,
								&C,
								STBI_rgb );
	std::vector<double> image_double(W*H*3);
	for (int i=0; i<W*H*3; i++)
		image_double[i] = image[i];

	for( int contor = 0; contor < ITER; contor++ ) {
		Vector direction = get_random_direction();
		std::vector<std::pair<double, int>> vec1( W * H );
		std::vector<std::pair<double, int>> vec2( W * H );

		for( int i = 0; i < H; i++ ) {
			for( int j = 0; j < W; j++ ) {
				int index = ( i * W + j ) * 3;
				Vector pixel1( image_double[index], image_double[index + 1], image_double[index + 2] );
				Vector pixel2( image2[index], image2[index + 1], image2[index + 2] );
				vec1[index / 3] = std::make_pair( dot( direction, pixel1 ), index / 3 );
				vec2[index / 3] = std::make_pair( dot( direction, pixel2 ), index / 3 );
			}
		}

		std::sort( vec1.begin(), vec1.end() );
		std::sort( vec1.begin(), vec2.end() );

		for( int i = 0; i < W * H; i ++ ) {
			int index1 = vec1[i].second;
			int index2 = vec2[i].second;

			Vector D = ( vec2[i].first - vec1[i].first ) * direction;
			image_double[ index1 * 3 + 0 ] += D[0];
			image_double[ index1 * 3 + 1 ] += D[1];
			image_double[ index1 * 3 + 2 ] += D[2];
		}

	}
	
	std::vector<unsigned char> image_result(W*H * 3, 0);
	// for (int i = 0; i < H; i++) {
	// 	for (int j = 0; j < W; j++) {

	// 		image_result[(i*W + j) * 3 + 0] = image_double[(i*W+j)*3+0]*0.5;
	// 		image_result[(i*W + j) * 3 + 1] = image_double[(i*W+j)*3+1]*0.3;
	// 		image_result[(i*W + j) * 3 + 2] = image_double[(i*W+j)*3+2]*0.2;
	// 	}
	// }
	for( int i = 0; i < W * H * 3; i++ ) {
		image_result[i] = std::min( 255, std::max( 0, image_double[i] ) );
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);
	stbi_image_free( image );
	stbi_image_free( image2 );
	std::cout<<"Success!"<<std::endl;

	return 0;
}