// #define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

template <class T>
inline T MIN(const T &a, const T &b) 
{  return a < b ? a : b; }

template <class T>
inline T MAX(const T &a, const T &b) 
{  return a > b ? a : b; }

inline int INDEX(int ix,int iy,int iz,int nx,int ny) 
{  return ix+nx*(iy+ny*iz); }

#define CELLO_STRING_LENGTH 255

