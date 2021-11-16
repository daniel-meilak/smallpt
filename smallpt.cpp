#include <cmath>   
#include <cstdlib> 
#include <cstdio>
#include <iostream>
#include <random>
#include <vector>

// // seed for randum num generator
// static unsigned long s[2] = {1,1};

// static inline unsigned long rotl(const unsigned long x, int k) {
// 	return (x << k) | (x >> (64 - k));
// }

// // random number generator, using xoroshiro128+ for speed (outputs in range 0-1)
// unsigned long next(void) {
// 	const unsigned long s0 = s[0];
// 	unsigned long s1 = s[1];
// 	const unsigned long result = s0 + s1;

// 	s1 ^= s0;
// 	s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
// 	s[1] = rotl(s1, 37); // c

// 	return (result >> 11) * 0x1.0p-53;
// }

unsigned short Xi[3] = {0,0,5};
double next(){
   return erand48(Xi);
}


// material types, used in radiance()
enum Refl_t{ DIFF, SPEC, REFR };

struct Vec{

   // position, also color (r,g,b)
   double x, y, z;

   // Constructor          
   Vec(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z){};

   // vector addition and subtraction
   Vec operator+(const Vec& b) const{ return Vec(x + b.x, y + b.y, z + b.z); }
   Vec operator-(const Vec& b) const{ return Vec(x - b.x, y - b.y, z - b.z); }

   // scalar multiplication
   Vec operator*(double b) const{ return Vec(x * b, y * b, z * b); }

   // vector multiplication
   Vec mult(const Vec& b) const{ return Vec(x * b.x, y * b.y, z * b.z); }

   // normalise vector
   Vec& norm(){ return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }

   // dot product
   double dot(const Vec& b) const{ return x * b.x + y * b.y + z * b.z; } // cross:

   // cross product
   Vec operator%(Vec& b){ return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Ray{
   Vec origin, direction;
   Ray(Vec origin, Vec direction): origin(origin), direction(direction){};
};

struct Sphere{

   double radius;
   Vec position, emission, color;

   // reflection type (DIFFuse, SPECular, REFRactive)
   Refl_t refl;

   // constructor
   Sphere(double radius, Vec position, Vec emission, Vec color, Refl_t refl):
      radius(radius), position(position), emission(emission), color(color), refl(refl){};

   // returns distance, 0 if nohit
   // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
   double intersect(const Ray& ray) const{
      Vec op = position - ray.origin;
      
      // fudge factor
      double t, eps = 1e-4;

      // 1/2 b from quadratic equation setup
      double b = op.dot(ray.direction);
      
      // (b^2 -4ac)/4 - a=1 due to normalization
      double det = b * b - op.dot(op) + radius * radius;

      // if ray misses sphere
      if (det < 0){ return 0; }
      else { det = sqrt(det); }

      return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
   }
};

std::vector<Sphere> spheres = {

   //   radius  position                  emission       color              material
   Sphere(1e5 , Vec(1e5 + 1,40.8,81.6)  , Vec()        , Vec(.75,.25,.25) , DIFF),  //Left
   Sphere(1e5 , Vec(-1e5 + 99,40.8,81.6), Vec()        , Vec(.25,.25,.75) , DIFF),  //Rght
   Sphere(1e5 , Vec(50,40.8, 1e5)       , Vec()        , Vec(.75,.75,.75) , DIFF),  //Back
   Sphere(1e5 , Vec(50,40.8,-1e5 + 170) , Vec()        , Vec()            , DIFF),  //Frnt
   Sphere(1e5 , Vec(50, 1e5, 81.6)      , Vec()        , Vec(.75,.75,.75) , DIFF),  //Botm
   Sphere(1e5 , Vec(50,-1e5 + 81.6,81.6), Vec()        , Vec(.75,.75,.75) , DIFF),  //Top
   Sphere(16.5, Vec(27,16.5,47)         , Vec()        , Vec(1,1,1) * .999, SPEC),  //Mirr
   Sphere(16.5, Vec(73,16.5,78)         , Vec()        , Vec(1,1,1) * .999, REFR),  //Glas
   Sphere(600 , Vec(50,681.6 - .27,81.6), Vec(12,12,12), Vec()            , DIFF)   //Lite
};

// value must be between 0-1
inline double clamp(double x){ return x < 0 ? 0 : x>1 ? 1 : x; }

// rgb value made by scaling to 0-255, with gamma correction of 2.2
inline int toInt(double x){ return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

// intersection of ray with scene
inline bool intersect(const Ray& ray, double& t, int& id){
   
   double d, inf = t = 1e20;
   
   for (size_t i = spheres.size(); i--;){
      
      if ((d = spheres[i].intersect(ray)) && d < t){
         t  = d;
         id = i;
      }
   }
   
   return t < inf;
}

Vec radiance(const Ray& ray, int depth){

   // distance to intersection
   double t;                  

   // id of intersected object
   int id = 0;

   // if miss, return black
   if (!intersect(ray, t, id)){ return Vec(); } 
   
   // the hit object
   const Sphere& obj = spheres[id];        
   
   // ray intersection point
   Vec x = ray.origin + ray.direction * t;

   // sphere normal
   Vec n = (x - obj.position).norm();
   
   // properly oriented surface normal
   Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;
   
   // object color (BRDF modulator)
   Vec f = obj.color;
   
   // use maximum reflectivity amount for russian roulette
   double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

   // russian roulette
   if (++depth > 5){
      if (next() < p){ f = f * (1 / p); }
      else { return obj.emission; }
   }
   
   // Ideal DIFFUSE reflection
   if (obj.refl == DIFF){                  
      
      // angle around
      double r1 = 2 * M_PI * next();
      
      // distance from center
      double r2 = next(), r2s = sqrt(r2);
      
      // normal
      Vec w = nl;

      // u is perpendicular to w
      Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
      
      // v is perpendicular to u and w
      Vec v = w % u;

      // d is random reflection ray
      Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
      
      return obj.emission + f.mult(radiance(Ray(x, d), depth));
   }
   // ideal SPECULAR reflection
   else if (obj.refl == SPEC){           
      return obj.emission + f.mult(radiance(Ray(x, ray.direction - n * 2 * n.dot(ray.direction)), depth));
   }

   // otherwise we have a dielectric (glass) surface -> ideal dielectic refraction
   Ray reflRay(x, ray.direction - n * 2 * n.dot(ray.direction));

   // ray from outside going in
   bool into = n.dot(nl) > 0; 

   double nc = 1;
   double nt = 1.5;
   double nnt = into ? nc / nt : nt / nc;
   double ddn = ray.direction.dot(nl);
   double cos2t;
   
   // Total internal reflection
   if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0){  
      return obj.emission + f.mult(radiance(reflRay, depth));
   }

   // otherwise, choose reflection or refraction
   Vec tdir = (ray.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
   double a = nt - nc, b = nt + nc;
   double R0 = a * a / (b * b);
   double c = 1 - (into ? -ddn : tdir.dot(n));
   double Re = R0 + (1 - R0) * c * c * c * c * c;
   double Tr = 1 - Re;
   double P = .25 + .5 * Re;
   double RP = Re / P;
   double TP = Tr / (1 - P);
   
   // Russian roulette
   return obj.emission + f.mult(depth > 2 ? (next() < P ?   
      radiance(reflRay, depth) * RP : radiance(Ray(x, tdir), depth) * TP):
      radiance(reflRay, depth) * Re + radiance(Ray(x, tdir), depth) * Tr);
}

int main(int argc, char* argv[]){
   
   // image size
   int width = 1024, height = 768;
   
   // samples (defaults to 1)
   int samps = argc == 2 ? atoi(argv[1]) / 4 : 1;
   
   // camera position and direction
   Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());

   // x direction increment
   Vec cx = Vec(width * .5135 / height);
   
   // y direction increment
   Vec cy = (cx % camera.direction).norm() * .5135;

   // Color of samples 
   Vec r;
   
   // image
   Vec* c = new Vec[width * height];
   
   // OpenMP
   #pragma omp parallel for schedule(dynamic, 1) private(r)       
   // Loop over image rows
   for (int y = 0; y < height; y++){

      // print progress
      // std::cerr << "\rRendering (" << samps * 4 << " spp) %" << 100. * y / (height - 1); 
      fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (height - 1));
      
      // Loop over columns
      for (int x = 0; x < width; x++)   
         
         // 2x2 subpixel rows
         for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++)
            
            // 2x2 subpixel cols
            for (int sx = 0; sx < 2; sx++, r = Vec()){        
               
               for (int s = 0; s < samps; s++){
                 
                  // Tent filter
                  double r1 = 2 * next(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                  double r2 = 2 * next(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                  
                  Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) + cy * (((sy + .5 + dy) / 2 + y) / height - .5) + camera.direction;
                  
                  r = r + radiance(Ray(camera.origin + d * 140, d.norm()), 0) * (1. / samps);
               } // Camera rays are pushed ^^^^^ forward to start in interior
               
               // Subpixel estimate
               c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
            }
   }
   std::cout << std::endl;

   // write image to ppm file
   FILE* f = fopen("image.ppm", "w");         
   fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
   for (int i = 0; i < width * height; i++)
      fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
