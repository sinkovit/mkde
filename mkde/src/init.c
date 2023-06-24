#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP computeCellSurfaceArea(void *, void *);
extern SEXP mkde2dGrid02(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP mkde2dGridv02interact(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP mkde3dGridv02(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP mkde3dGridv02interact(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP writeMKDE3DtoGRASS(void *, void *, void *, void *, void *, void *);
extern SEXP writeMKDE3DtoVTK(void *, void *, void *, void *, void *, void *);
extern SEXP writeMKDE3DtoXDMF(void *, void *, void *, void *, void *, void *);
extern SEXP writeRasterToVTK02(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP writeRasterToXDMF02(void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"computeCellSurfaceArea", (DL_FUNC) &computeCellSurfaceArea,  2},
    {"mkde2dGrid02",           (DL_FUNC) &mkde2dGrid02,           10},
    {"mkde2dGridv02interact",  (DL_FUNC) &mkde2dGridv02interact,  14},
    {"mkde3dGridv02",          (DL_FUNC) &mkde3dGridv02,          16},
    {"mkde3dGridv02interact",  (DL_FUNC) &mkde3dGridv02interact,  23},
    {"writeMKDE3DtoGRASS",     (DL_FUNC) &writeMKDE3DtoGRASS,      6},
    {"writeMKDE3DtoVTK",       (DL_FUNC) &writeMKDE3DtoVTK,        6},
    {"writeMKDE3DtoXDMF",      (DL_FUNC) &writeMKDE3DtoXDMF,       6},
    {"writeRasterToVTK02",     (DL_FUNC) &writeRasterToVTK02,    8},
    {"writeRasterToXDMF02",    (DL_FUNC) &writeRasterToXDMF02,   5},
    {NULL, NULL, 0}
};

void R_init_mkde(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
