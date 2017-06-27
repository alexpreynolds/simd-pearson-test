#ifndef PEARSON_H_
#define PEARSON_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif /* getline() support */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdbool.h>
#include <immintrin.h>
#include <errno.h>

#define ID_MAX_LEN 524288
#define FN_MAX_LEN 1024
#define BUF_MAX_LEN 4096
#define CHR_MAX_LEN 256
#define COORD_MAX_LEN 20
#define SCORE_MAX_LEN 5
#define DELIMITER_LEN 1
#define ENTRY_MAX_LEN 20
#define AVX_FLOAT_N 8
#define ALIGNMENT 32

    typedef unsigned char byte_t;
    typedef float score_t;

    extern const int kSignalDelim;
    const int kSignalDelim = (int) ',';

    extern const byte_t kNANEncodedByte;
    const byte_t kNANEncodedByte = 0xca;
    
    extern const score_t kEpsilon;
    const score_t kEpsilon = 0.0000001f;

    extern const score_t kEpsilonLessStringent;
    const score_t kEpsilonLessStringent = 0.001f;
    
    extern const score_t kSelfCorrelationScore;
    const score_t kSelfCorrelationScore = +1.0f;
    
    extern const score_t kNoCorrelationScore;
    const score_t kNoCorrelationScore = +0.0f;

    static const score_t pt_encode_byte_to_score_table[256] = 
        {-1.00, 
         -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.90,
         -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.80,
         -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.70, 
         -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.60, 
         -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.50, 
         -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.40, 
         -0.39, -0.38, -0.37, -0.36, -0.35, -0.34, -0.33, -0.32, -0.31, -0.30, 
         -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23, -0.22, -0.21, -0.20, 
         -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, 
         -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, -0.00,
         +0.00, +0.01, +0.02, +0.03, +0.04, +0.05, +0.06, +0.07, +0.08, +0.09,
         +0.10, +0.11, +0.12, +0.13, +0.14, +0.15, +0.16, +0.17, +0.18, +0.19,
         +0.20, +0.21, +0.22, +0.23, +0.24, +0.25, +0.26, +0.27, +0.28, +0.29,
         +0.30, +0.31, +0.32, +0.33, +0.34, +0.35, +0.36, +0.37, +0.38, +0.39,
         +0.40, +0.41, +0.42, +0.43, +0.44, +0.45, +0.46, +0.47, +0.48, +0.49,
         +0.50, +0.51, +0.52, +0.53, +0.54, +0.55, +0.56, +0.57, +0.58, +0.59,
         +0.60, +0.61, +0.62, +0.63, +0.64, +0.65, +0.66, +0.67, +0.68, +0.69,
         +0.70, +0.71, +0.72, +0.73, +0.74, +0.75, +0.76, +0.77, +0.78, +0.79,
         +0.80, +0.81, +0.82, +0.83, +0.84, +0.85, +0.86, +0.87, +0.88, +0.89,
         +0.90, +0.91, +0.92, +0.93, +0.94, +0.95, +0.96, +0.97, +0.98, +0.99,
         +1.00, 
           NAN, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, 
         +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, 
         +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, 
         +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, 
         +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, 
         +0.00, +0.00, +0.00, +0.00};

    typedef struct signal {
        uint32_t n;
        score_t* data;
        score_t mean;
        score_t sd;
    } signal_t;

    typedef struct bed5 {
        char* chr;
        uint64_t start;
        uint64_t stop;
        char* id;
        signal_t* signal;
    } bed5_t;

    typedef struct signal_avx {
        uint32_t n_raw;
        uint32_t n;
        __m256* data;
        score_t mean;
        score_t sd;
    } signal_avx_t;

    typedef struct bed5_avx {
        char* chr;
        uint64_t start;
        uint64_t stop;
        char* id;
        signal_avx_t* signal;
    } bed5_avx_t;

    typedef struct lookup {
        uint64_t capacity;
        uint32_t nelems;
        bed5_t** elems;
    } lookup_t;

    typedef struct lookup_avx {
        uint64_t capacity;
        uint32_t nelems;
        bed5_avx_t** elems;
    } lookup_avx_t;

    static struct pt_globals_t {
        bool verbose;
        char* input_fn;
        char* output_fn;
        FILE* output_stream;
        bool avx_enabled;
        lookup_t* lookup;
        lookup_avx_t* lookup_avx;
    } pt_globals;

    void pt_initialize_lookup_avx_via_signal_avx(const char* ifn, lookup_avx_t** lp);    
    void pt_initialize_lookup_via_signal(const char* ifn, lookup_t** lp);
    void pt_initialize_bed5_avx_element(char* chr, uint64_t start, uint64_t stop, char* id, bed5_avx_t** ep);
    void pt_initialize_bed5_element(char* chr, uint64_t start, uint64_t stop, char* id, bed5_t** ep);
    void pt_initialize_signal_avx(char* id, signal_avx_t** sp);
    void pt_initialize_signal(char* id, signal_t** sp);
    void pt_push_bed5_avx_element_to_lookup_avx(bed5_avx_t* e, lookup_avx_t** lp);
    void pt_push_bed5_element_to_lookup(bed5_t* e, lookup_t** lp);
    void pt_delete_lookup_avx(lookup_avx_t** lp);
    void pt_delete_lookup(lookup_t** lp);
    void pt_delete_bed5_avx_element(bed5_avx_t** ep);
    void pt_delete_bed5_element(bed5_t** ep);
    void pt_delete_signal_avx(signal_avx_t** sp);
    void pt_delete_signal(signal_t** sp);
    void pt_initialize_globals();
    void pt_initialize_command_line_options(int argc, char** argv);
    void pt_print_usage(FILE* os);
    void pt_delete_globals();
    score_t pt_pearson_r_via_signal_t(signal_t* a, signal_t* b);
    static inline score_t pt_mean_signal_avx(__m256* d, uint32_t len);
    score_t pt_mean_signal(score_t* d, uint32_t len);
    score_t pt_sample_sd_signal_avx(score_t* d, uint32_t len, score_t m);
    score_t pt_sample_sd_signal(score_t* d, uint32_t len, score_t m);
    inline score_t pt_truncate_score_to_precision(score_t d, int prec);
    inline byte_t pt_encode_score_to_byte(score_t d);
    inline bool pt_signbit(score_t d);
    void pt_calculate_pearson_scores_via_signal_avx(lookup_avx_t* l);
    void pt_calculate_pearson_scores_via_signal(lookup_t* l);

    static const char* pt_client_opt_string = "ai:o:vh?";

    static struct option pt_client_long_options[] = {
        { "avx",                                        no_argument,       NULL, 'a' },
        { "input",                                      required_argument, NULL, 'i' },
        { "output",                                     required_argument, NULL, 'o' },
        { "verbose",                                    no_argument,       NULL, 'v' },
        { "help",                                       no_argument,       NULL, 'h' },
        { NULL,                                         no_argument,       NULL,  0  }
    };

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif // PEARSON_H_
