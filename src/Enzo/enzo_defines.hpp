#ifndef ENZO_DEFINES_H
#define ENZO_DEFINES_H

#define ENZO_FAIL     0
#define ENZO_SUCCESS  1

#define FALSE         0
#define TRUE          1

#define ENZO_FLOAT_UNDEFINED  -99999.0
#define ISYM         "d"

#define ENZO_INDEX_UNDEFINED -1

#define MAX_DIMENSION                3
#define MAX_NUMBER_OF_BARYON_FIELDS 100

#define MAX_COLOR                20
#define COLOR_FLOOR              1.d-35
#define MAX_ANY_SINGLE_DIRECTION (1024+6)

#define NUM_FIELDS 34

#define NINT(A)   ((int) ((A) + 0.5*((A)>0?1:-1)))

#define FORTRAN_NAME(NAME) NAME##_

#define Density          0
#define TotalEnergy      1
#define InternalEnergy   2
#define Velocity1        4
#define Velocity2        5
#define Velocity3        6
#define ElectronDensity  7
#define HIDensity        8
#define HIIDensity       9
#define HeIDensity       10
#define HeIIDensity      11
#define HeIIIDensity     12
#define HMDensity        13
#define H2IDensity       14
#define H2IIDensity      15
#define DIDensity        16
#define DIIDensity       17
#define HDIDensity       18
#define Metallicity      19

#endif /* ENZO_DEFINES_HPP */
