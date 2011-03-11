// #define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

template <class T>
T MIN(const T &a, const T &b) {  return a < b ? a : b; }

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define INDEX(ix,iy,iz,nx,ny) ((ix)+(nx)*((iy)+(ny)*(iz)))

