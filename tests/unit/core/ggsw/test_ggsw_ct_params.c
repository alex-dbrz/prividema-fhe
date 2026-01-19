#include "core/ggsw/ggsw_ct_params.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE 10
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)

Test(new_ggsw_ct_params, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, params->k_tilde, K_TILDEBASE));
    cr_assert(eq(i64, params->kappa_tilde, KAPPA_TILDEBASE));
    cr_assert(eq(i64, params->n_limbs_tilde, NLIMBS_TILDEBASE));

    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}