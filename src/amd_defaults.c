/* ========================================================================= */
/* === AMD_defaults ======================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.                                     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable.  Sets default control parameters for AMD.  See amd.h
 * for details.
 */

#include "amd_internal.h"

/* ========================================================================= */
/* === AMD defaults ======================================================== */
/* ========================================================================= */

GLOBAL void AMD_defaults
(
    double Control [ ]
)
{
    Int i ;

    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < AMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [AMD_DENSE] = AMD_DEFAULT_DENSE ;
	Control [AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE ;
    }
}
