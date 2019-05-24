#ifndef	__LES_LOCAL_H__
#define	__LES_LOCAL_H__
#include <FCMangle.h>
/*===========================================================================
 *
 * "Top":  A dummy top driver
 *
 *===========================================================================
 */
typedef struct _Top* TopHd ;
typedef struct _Usr* UsrHd ;

/*===========================================================================
 *
 * "Les":  Les data structure
 *
 *===========================================================================
 */
typedef struct _Les {
    TopHd	topHd ;			/* top driver			*/

    Real*	stats ;			/* solution statistics		*/

    Integer	nPresPrjs ;		/* No. pres prjs		*/
    Integer	nVelPrjs ;		/* No. vel  prjs		*/
    Integer	nFlowPrjs ;		/* No. flow prjs		*/
    Integer	nTempPrjs ;		/* No. temp prjs		*/
    Integer	nTurbPrjs ;		/* No. turb prjs		*/
    Integer	nSpecPrjs ;		/* No. spec prjs		*/

    Integer	nActPresPrjs ;		/* No. pres prjs current	*/
    Integer	nActVelPrjs ;		/* No. vel  prjs current	*/
    Integer	nActFlowPrjs ;		/* No. flow prjs current	*/
    Integer	nActTempPrjs ;		/* No. temp prjs current	*/
    Integer	nActTurbPrjs ;		/* No. turb prjs current	*/
    Integer*	nActSpecPrjs ;		/* No. spec prjs current	*/

    Real	presRegFct ;		/* pressure regularization	*/
    Real	velRegFct ;		/* velocity regularization	*/
    Real	tempRegFct ;		/* temperature regularization	*/
    Real	turbRegFct ;		/* turbulence regularization	*/
    Real	specRegFct ;		/* species regularization	*/

    Real*	presCoef ;		/* pres prj coefs		*/
    Real*	velCoef ;		/* vel  prj coefs		*/
    Real*	flowCoef ;		/* flow prj coefs		*/
    Real*	tempCoef ;		/* temp prj coefs		*/
    Real*	turbCoef ;		/* turb prj coefs		*/
    Real**	specCoef ;		/* spec prj coefs		*/

    Real*	memTmp ;		/* temporary memory		*/
    Real*	memTmpCurr ;		/* memTmp current position	*/
    Integer	memTmpDim ;		/* dimension of tmpMem		*/

    Integer	nTmpVecs ;		/* No. temporary vectors	*/
    Integer	nPermVecs ;		/* No. Permanent vectors	*/

    Integer	mPresS ;		/* pres S vector		*/
    Integer	mVelS ;			/* vel  S vector		*/
    Integer	mFlowS ;		/* flow S vector		*/
    Integer	mTempS ;		/* temp S vector		*/
    Integer	mTurbS ;		/* turb S vector		*/
    Integer*	mSpecS ;		/* spec S vector		*/
    Integer	mRes ;			/* res vector			*/
    Integer	mSol ;			/* solInc vector		*/

    Integer	eqnType ;		/* equation being solved	*/
    Integer	specId ;		/* species being solved		*/
    Integer	lesType ;		/* linear solver		*/
    Real	tol ;			/* convergence tolerance	*/
    Real	prjTol ;		/* convergence tolerance	*/
    Integer	minIters ;		/* Min No. iterations		*/
    Integer	maxIters ;		/* Max No. iterations		*/
    Integer	nKvecs ;		/* No. Krylov Vectors		*/
    Integer	prjFlag ;		/* Projection flag		*/
    Integer	presPrjFlag ;		/* Pressure projection flag	*/
    Integer	presPrecFlag ;		/* Pressure precondition flag	*/
    Integer	mlpFlag ;		/* multi-level precond. flag	*/

    Integer	nDofs ;			/* No. dofs			*/

    Integer	verbose ;		/* verbose level		*/
} Les ;

typedef struct _Les* LesHd ;
/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */

#endif	/* __LES_LOCAL_H__ */
