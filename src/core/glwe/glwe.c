#include "glwe.h"
#include "vec_znx_arithmetic_private.h"
#include "distributions.h"

/**
 * @brief Encrypts the phase (message + noise) and puts it in res.
 * 
 * @param module The module stocking the degree N.
 * @param params The GLWE parameters.
 * @param res_ct The result bivariate ciphertext. 
 * @param sk The secret key in DFT space.
 * @param phase message + noise.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int glwe_encrypt_priv(const Core* core,    // all params of the library: is fft
                                            // or ntt, all N that are used
                       GLWECiphertext* ct,  // ciphertext
                       GLWEPreparedSK* sk_dft,   // secret key: vec of size k
                       PolyBiv* phase    // message + noise
){
    int64_t N = ct->params->N;
    int64_t k = ct->params->k;
    int64_t kappa = ct->params->kappa;
    int64_t n_limbs = ct->params->n_limbs;
    int64_t l = n_limbs / (k+1);

    MODULE* module = new_module_info(N, FFT64)
    if (uniform_random_vec(k * N, ct->vec, l, (k + 1) * N, kappa) > 0) {
        return -1;
    }
    
    // acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes ∑_j{0,k-1}[s_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes DFT(s_j) * DFT(a_j)
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, resVec_j_dft, l, sk_j_univ_dft, ct->vec + j*N, l, (k+1)*N); 
        
        // Computes s_j * a_j
        PolyBiv* resVec_j = new_vec_znx_big_p(module, l); 
        vec_znx_idft_p(module, resVec_j, l, resVec_j_dft, l);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] += resVec_j[p];
        }
        delete_vec_znx_dft_p(resVec_j_dft);
        delete_vec_znx_big_p(resVec_j);
    }

    // Add the phase to acc
    for(int64_t i = 0 ; i < l ; i++){    
        for(int64_t p = 0 ; p < N ; p++){
                acc[i*N + p] += phase->coeffs[i*N + p];
        }
    }
    
    // The pointer to limb_0(b)
    PolyBiv* b_0 = ct->vec + k*N;

    // For each i in {0,l} limb_i(b) = limb_i(acc) = ∑_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k_p(module, kappa, b_0, l, N*(k+1), acc, l, N);
    
    free(acc);

    return 0;
}

/**
 * @brief Decrypts the phase (message + noise) and puts it in phase.
 * 
 * @param enc_params The GLWE parameters.
 * @param phase The phase in Rn[X]. 
 * @param key The secret key in DFT space.
 * @param ct The ciphertext.
 */
int decrypt_biv_glwe(GLWECtParams* enc_params,
                     TNXElement* phase, 
                     GLWEPreparedSK* sk_dft,
                     GLWECiphertext* ct
){
    // GLWE parameters
    int64_t N = enc_params->N;
    int64_t k = enc_params->k;
    int64_t l = poly_biv_size(enc_params);

    MODULE* module = new_module_info(N, FFT64);

    PolyBiv* acc = malloc(poly_biv_bytes(enc_params)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes acc = -∑_j{0,k-1}[sk_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of resp. the secret key and the bivGLWE ciphertext 
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        PolyUniv* a_j = ct->vec + j*N;
        
        // Computes DFT(sk_j * a_j)
        PolyBivDFT* as_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k+1)*N); 
        
        // Computes sk_j * a_j
        PolyBiv* as_j = new_vec_znx_big_p(module, l); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // And subs it to acc
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] -= as_j[p];
        }
        delete_vec_znx_dft_p(as_j_dft);
        delete_vec_znx_big_p(as_j);
    }

    // Computes acc = b - ∑_j{0,k-1}[sk_j * a_j]
    int64_t* b = ct->vec + N*k;
    add_biv_poly(enc_params, acc, N, b, N*(k+1), acc, N);
    
    biv_to_univ(enc_params, phase, acc);

    return 0;
}
