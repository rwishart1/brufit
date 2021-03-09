 #ifndef PLANCK_ALLOC_UTILS_H
 #define PLANCK_ALLOC_UTILS_H
 
 #include <cstdlib>
 #include <cstddef>
 #include "datatypes.h"
 
 template <typename T> class normalAlloc__
   {
   public:
     static T *alloc(tsize sz) { return (sz>0) ? new T[sz] : 0; }
     static void dealloc (T *ptr) { delete[] ptr; }
   };
 
 template <typename T, int align> class alignAlloc__
   {
   private:
 //#if (__cplusplus>=201103L)
 //    enum { max_nat_align = alignof(std::max_align_t) };
 //#else
     enum { max_nat_align = sizeof(void *) };
 //#endif
 
   public:
     static T *alloc(tsize sz)
       {
       using namespace std;
       if (sz==0) return 0;
       planck_assert((align&(align-1))==0,"alignment must be power of 2");
       void *res;
       if (align<=max_nat_align)
         {
         res=malloc(sz*sizeof(T));
         planck_assert(res!=0,"error in malloc()");
         }
       else
         {
         tsize overhead=align-1+sizeof(void*);
         void *ptr=malloc(sz*sizeof(T)+overhead);
         planck_assert(ptr!=0,"error in malloc()");
         tsize sptr=reinterpret_cast<tsize>(ptr);
         sptr = (sptr+overhead) & ~(align-1);
         void **ptr2 = reinterpret_cast<void **>(sptr);
         ptr2[-1]=ptr;
         res=ptr2;
         }
       return static_cast<T *>(res);
       }
     static void dealloc(T *ptr)
       {
       using namespace std;
       if (align<=max_nat_align)
         free(ptr);
       else
         {
         if (ptr==0) return;
         void **ptr2 = reinterpret_cast<void **>(ptr);
         free (ptr2[-1]);
         }
       }
   };
 
 #endif
