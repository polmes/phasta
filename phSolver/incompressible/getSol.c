#include "usr.h"
/*#include "rts.h" */
void getSol ( UsrHd usrHd,
#ifdef SP_Solve
              float* Dy  )
#else
              double* Dy  )
#endif
{

     Dy = usrHd->solinc;
     
/* extern int rts_getavailablememory(size_t* availableMemory);
*/

}
