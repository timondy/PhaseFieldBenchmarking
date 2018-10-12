#include <string.h>
#include "underscore.h"

#ifdef UNDERSCORE
void c_get_pointer_ (long b, long *ptr) {
#endif
#ifdef DOUBLE_UNDERSCORE
void c_get_pointer__ (long b, long *ptr) {
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
void c_get_pointer (long b, long *ptr) {
#endif
#endif

  *ptr = b;

}

#ifdef UNDERSCORE
extern void fortran_get_pointer_();
#endif
#ifdef DOUBLE_UNDERSCORE
extern void fortran_get_pointer__();
#endif
#ifndef UNDERSCORE
#ifndef DOUBLE_UNDERSCORE
extern void fortran_get_pointer();
#endif
#endif

void* c_amr_get_pointer(char* name) {

   void* ptr;
   int istrlen;

   istrlen = strlen(name);

   ptr = NULL;

#ifdef UNDERSCORE
   fortran_get_pointer_(&ptr,name,&istrlen);
#endif
#ifdef DOUBLE_UNDERSCORE
   fortran_get_pointer__(&ptr,name,&istrlen);
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
   fortran_get_pointer(&ptr,name,&istrlen);
#endif
#endif   

   return(ptr);
}
