// #define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

template <class T>
T MIN(const T &a, const T &b) {  return a < b ? a : b; }

template <class T>
T MAX(const T &a, const T &b) {  return a > b ? a : b; }

#define INDEX(ix,iy,iz,nx,ny) ((ix)+(nx)*((iy)+(ny)*(iz)))

#define CELLO_STRING_LENGTH 255

