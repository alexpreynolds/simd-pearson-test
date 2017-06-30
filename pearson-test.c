#include "pearson-test.h"

int
main(int argc, char** argv) 
{
    bs_initialize_globals();
    bs_initialize_command_line_options(argc, argv);
    if (!bs_globals.avx_enabled) {
        bs_initialize_lookup_via_signal(bs_globals.input_fn, &bs_globals.lookup);
        bs_calculate_pearson_scores_via_signal(bs_globals.lookup);
    }
    else {
        bs_initialize_lookup_avx_via_signal_avx(bs_globals.input_fn, &bs_globals.lookup_avx);
        bs_calculate_pearson_scores_via_signal_avx(bs_globals.lookup_avx);
    }
    bs_delete_globals();

    return EXIT_SUCCESS;
}

void
bs_initialize_lookup_avx_via_signal_avx(const char* ifn, lookup_avx_t** lp)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_lookup_avx_via_signal_avx()\n");

    lookup_avx_t* l = NULL;
    FILE* lf = NULL;
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char chr_str[CHR_MAX_LEN] = {0};
    char start_str[COORD_MAX_LEN] = {0};
    char stop_str[COORD_MAX_LEN] = {0};
    char id_str[ID_MAX_LEN] = {0};
    uint64_t start_val = 0;
    uint64_t stop_val = 0;
    
    l = malloc(sizeof(lookup_avx_t));
    if (!l) {
        fprintf(stderr, "Error: Could not allocate space for lookup table!\n");
        exit(EXIT_FAILURE);
    }
    l->capacity = 0;
    l->nelems = 0;
    l->elems = NULL;

    lf = fopen(ifn, "r");
    if (ferror(lf)) {
        fprintf(stderr, "Error: Could not open or read from [%s]\n", ifn);
        exit(EXIT_FAILURE);
    }

    /* parse BED element into bed5_t* and push to lookup table */
    while ((buf_read = getline(&buf, &buf_len, lf)) != -1) {
        sscanf(buf, "%s\t%s\t%s\t%s\n", chr_str, start_str, stop_str, id_str);
        sscanf(start_str, "%" SCNu64, &start_val);
        sscanf(stop_str, "%" SCNu64, &stop_val);
        bed5_avx_t* e = NULL;
        signal_avx_t* sa = NULL;
        //fprintf(stderr, "processing buf [%s] [%u] [%u] [%s]\n", chr_str, (uint32_t) start_val, (uint32_t) stop_val, id_str);
        bs_initialize_bed5_avx_element(chr_str, start_val, stop_val, id_str, sa, &e);
        bs_push_bed5_avx_element_to_lookup_avx(e, &l);
    }

    free(buf);
    fclose(lf);

    if (!l) {
        fprintf(stderr, "Error: Could not parse BED elements into lookup_avx struct\n");
        exit(EXIT_FAILURE);
    }

    *lp = l;
}

void
bs_initialize_lookup_via_signal(const char* ifn, lookup_t** lp)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_lookup_via_signal()\n");

    lookup_t* l = NULL;
    FILE* lf = NULL;
    char* buf = NULL;
    size_t buf_len = 0;
    ssize_t buf_read = 0;
    char chr_str[CHR_MAX_LEN] = {0};
    char start_str[COORD_MAX_LEN] = {0};
    char stop_str[COORD_MAX_LEN] = {0};
    char id_str[ID_MAX_LEN] = {0};
    uint64_t start_val = 0;
    uint64_t stop_val = 0;
    
    l = malloc(sizeof(lookup_t));
    if (!l) {
        fprintf(stderr, "Error: Could not allocate space for lookup table!\n");
        exit(EXIT_FAILURE);
    }
    l->capacity = 0;
    l->nelems = 0;
    l->elems = NULL;

    lf = fopen(ifn, "r");
    if (ferror(lf)) {
        fprintf(stderr, "Error: Could not open or read from [%s]\n", ifn);
        exit(EXIT_FAILURE);
    }

    /* parse BED element into bed5_t* and push to lookup table */
    while ((buf_read = getline(&buf, &buf_len, lf)) != -1) {
        sscanf(buf, "%s\t%s\t%s\t%s\n", chr_str, start_str, stop_str, id_str);
        sscanf(start_str, "%" SCNu64, &start_val);
        sscanf(stop_str, "%" SCNu64, &stop_val);
        bed5_t* e = NULL;
        bs_initialize_bed5_element(chr_str, start_val, stop_val, id_str, &e);
        bs_push_bed5_element_to_lookup(e, &l);
    }

    free(buf);
    fclose(lf);

    if (!l) {
        fprintf(stderr, "Error: Could not parse BED elements into lookup struct\n");
        exit(EXIT_FAILURE);
    }

    *lp = l;
}

void
bs_initialize_bed5_avx_element(char* chr, uint64_t start, uint64_t stop, char* id, signal_avx_t* sa, bed5_avx_t** ep)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_bed5_avx_element() [%p]\n", sa);

    bed5_avx_t *e = NULL;
    e = malloc(sizeof(bed5_avx_t));
    if (!e) {
        fprintf(stderr, "Error: Could not allocate space for bed5_avx element!\n");
        exit(EXIT_FAILURE);
    }
    else {
        //fprintf(stderr, "allocated new element\n");
    }
    
    e->chr = NULL;
    if (strlen(chr) > 0) {
        e->chr = malloc(strlen(chr) + 1);
        if (!e->chr) {
            fprintf(stderr, "Error: Could not allocate space for bed5_avx element chromosome!\n");
            exit(EXIT_FAILURE);
        }
        memcpy(e->chr, chr, strlen(chr) + 1);
        //fprintf(stderr, "copied chromosome name to new element\n");
    }
    
    e->start = start;
    //fprintf(stderr, "copied start position to new element\n");
    
    e->stop = stop;
    //fprintf(stderr, "copied stop position to new element\n");
    
    e->id = NULL;
    if (strlen(id) > 0) {
        e->id = malloc(sizeof(*id) * strlen(id) + 1);
        if (!e->id) {
            fprintf(stderr,"Error: Could not allocate space for bed5_avx element id!\n");
            exit(EXIT_FAILURE);
        }
        memcpy(e->id, id, strlen(id) + 1);
        //fprintf(stderr, "copied ID to new element\n");
    }
    else {
        fprintf(stderr, "could not copy ID to new element\n");
    }
    
    e->signal = NULL;
    if (!sa) { 
        //fprintf(stderr, "setting up new signal()\n");
        bs_initialize_signal_avx(e->id, &(e->signal));
    }
    else {
        //fprintf(stderr, "copying old signal()\n");
        bs_copy_signal_avx(sa, &(e->signal));
    }

    if (!e) {
        fprintf(stderr, "Error: Could not parse BED elements into bed5_avx_t struct\n");
        exit(EXIT_FAILURE);
    }
    
    *ep = e;
}

void
bs_initialize_bed5_element(char* chr, uint64_t start, uint64_t stop, char* id, bed5_t** ep)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_bed5_element()\n");

    bed5_t *e = NULL;
    e = malloc(sizeof(bed5_t));
    if (!e) {
        fprintf(stderr, "Error: Could not allocate space for bed5 element!\n");
        exit(EXIT_FAILURE);
    }
    e->chr = NULL;
    if (strlen(chr) > 0) {
        e->chr = malloc(sizeof(*chr) * strlen(chr) + 1);
        if (!e->chr) {
            fprintf(stderr, "Error: Could not allocate space for bed5 element chromosome!\n");
            exit(EXIT_FAILURE);
        }
        memcpy(e->chr, chr, strlen(chr) + 1);
    }
    e->start = start;
    e->stop = stop;
    e->id = NULL;
    if (id && (strlen(id) > 0)) {
        e->id = malloc(sizeof(*id) * strlen(id) + 1);
        if (!e->id) {
            fprintf(stderr,"Error: Could not allocate space for bed5 element id!\n");
            exit(EXIT_FAILURE);
        }
        memcpy(e->id, id, strlen(id) + 1);
    }
    (e->id) ? bs_initialize_signal(e->id, &(e->signal)) : NULL;

    if (!e) {
        fprintf(stderr, "Error: Could not parse BED elements into bed5_t struct\n");
        exit(EXIT_FAILURE);
    }
    
    *ep = e;
}

void
bs_copy_signal_avx(signal_avx_t* src, signal_avx_t** dest)
{
    signal_avx_t* s = NULL;
    s = malloc(sizeof(signal_avx_t));
    if (!s) {
        fprintf(stderr, "Error: Could not allocate space for signal_avx pointer copy!\n");
        exit(EXIT_FAILURE);
    }
    s->n = src->n;
    s->data = NULL;
    s->data = aligned_alloc(32, s->n * sizeof(*s->data));
    if (!s->data) {
        fprintf(stderr, "Error: Could not allocate space for signal_avx data pointer!\n");
        exit(EXIT_FAILURE);
    }
    for (uint32_t idx = 0; idx < s->n; idx++) {
        s->data[idx] = src->data[idx];
    }
    s->mean = src->mean;
    s->sd = src->sd;

    *dest = s;
}

void
bs_initialize_signal_avx(char* id, signal_avx_t** sp)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_signal_avx()\n");

    signal_avx_t* s = NULL;
    s = malloc(sizeof(signal_avx_t));
    if (!s) {
        fprintf(stderr, "Error: Could not allocate space for signal_avx pointer!\n");
        exit(EXIT_FAILURE);
    }

    s->n = 0;
    s->data = NULL;
    s->mean = NAN;
    s->sd = NAN;

    for (uint32_t idx = 0; idx < strlen(id); idx++) {
        if (id[idx] == kSignalDelim) {
            s->n++;
        }
    }
    s->n++; // one more element than number of delimiters

    s->data = aligned_alloc(32, s->n * sizeof(*s->data));
    if (!s->data) {
        fprintf(stderr, "Error: Could not allocate space for signal data pointer!\n");
        exit(EXIT_FAILURE);
    }

    char* start = id;
    char* end = id;
    char entry_buf[ENTRY_MAX_LEN];
    uint32_t entry_idx = 0;
    bool finished_parsing = false;
    bool data_contains_nan = false;
    do {
        end = strchr(start, kSignalDelim);
        if (!end) {
            end = id + strlen(id);
            finished_parsing = true;
        }
        memcpy(entry_buf, start, end - start);
        entry_buf[end - start] = '\0';
        sscanf(entry_buf, "%f", &s->data[entry_idx++]);
        if (isnan(s->data[entry_idx - 1])) {
            data_contains_nan = true;
        }
        start = end + 1;
    } while (!finished_parsing);

    if (!data_contains_nan && s->n >= 2) {
        bs_calculate_mean_and_sd_signal_avx(s->data, s->n, &s->mean, &s->sd);
    }
    else {
        fprintf(stderr, "Warning: Vector has one value and therefore cannot have a standard deviation!\n");
    }

    *sp = s;
}

void
bs_initialize_signal(char* id, signal_t** sp)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_initialize_signal()\n");

    signal_t* s = NULL;
    s = malloc(sizeof(signal_t));
    if (!s) {
        fprintf(stderr, "Error: Could not allocate space for signal pointer!\n");
        exit(EXIT_FAILURE);
    }
    s->n = 0;
    s->data = NULL;
    s->mean = NAN;
    s->sd = NAN;

    for (uint32_t idx = 0; idx < strlen(id); idx++) {
        if (id[idx] == kSignalDelim) {
            s->n++;
        }
    }
    s->n++;

    s->data = malloc(sizeof(*s->data) * s->n);
    if (!s->data) {
        fprintf(stderr, "Error: Could not allocate space for signal data pointer!\n");
        exit(EXIT_FAILURE);
    }
    char* start = id;
    char* end = id;
    char entry_buf[ENTRY_MAX_LEN];
    uint32_t entry_idx = 0;
    bool finished_parsing = false;
    bool data_contains_nan = false;
    do {
        end = strchr(start, kSignalDelim);
        if (!end) {
            end = id + strlen(id);
            finished_parsing = true;
        }
        memcpy(entry_buf, start, end - start);
        entry_buf[end - start] = '\0';
        sscanf(entry_buf, "%f", &s->data[entry_idx++]);
        if (isnan(s->data[entry_idx - 1])) {
            data_contains_nan = true;
        }
        start = end + 1;
    } while (!finished_parsing);

    if (!data_contains_nan) {
        s->mean = bs_mean_signal(s->data, s->n);
        if (s->n >= 2) {
            s->sd = bs_sample_sd_signal(s->data, s->n, s->mean);
        }
        else {
            fprintf(stderr, "Warning: Vector has one value and therefore cannot have a standard deviation!\n");
        }
    }

    *sp = s;
}

void
bs_push_bed5_avx_element_to_lookup_avx(bed5_avx_t* e, lookup_avx_t** lp)
{
    if (bs_globals.verbose) fprintf(stderr, "-----------\nbs_push_bed5_avx_element_to_lookup_avx()\n");

    if ((*lp)->capacity == 0) {
        (*lp)->capacity++;
        (*lp)->elems = malloc(sizeof(bed5_avx_t *));
    }
    else if ((*lp)->nelems >= (*lp)->capacity) {
        if (bs_globals.verbose) fprintf(stderr, "resizing lookup table [capacity %d -> %d]...\n", (int)(*lp)->capacity, (int)(*lp)->capacity * 2);
        (*lp)->capacity *= 2;
        bed5_avx_t** new_elems = malloc(sizeof(bed5_avx_t *) * (*lp)->capacity);
        for (uint32_t idx = 0; idx < (*lp)->nelems; idx++) {
            if (bs_globals.verbose) fprintf(stderr, "copying old element to new at index [%u]...\n", idx + 1);
            bs_initialize_bed5_avx_element((*lp)->elems[idx]->chr,
                                           (*lp)->elems[idx]->start,
                                           (*lp)->elems[idx]->stop,
                                           (*lp)->elems[idx]->id,
                                           (*lp)->elems[idx]->signal,
                                           &new_elems[idx]);
            if (bs_globals.verbose) fprintf(stderr, "deleting old element at index [%u]...\n", idx + 1);
            bs_delete_bed5_avx_element(&((*lp)->elems[idx]));
        }   
        free((*lp)->elems);
        (*lp)->elems = NULL;
        (*lp)->elems = new_elems;     
    }
    if (bs_globals.verbose) fprintf(stderr, "number of elements now [%d]\n", (*lp)->nelems + 1);
    uint32_t n = (*lp)->nelems;
    (*lp)->elems[n] = e;
    (*lp)->nelems++;
}

void
bs_push_bed5_element_to_lookup(bed5_t* e, lookup_t** lp)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_push_bed5_element_to_lookup()\n");

    if ((*lp)->capacity == 0) {
        (*lp)->capacity++;
        (*lp)->elems = malloc(sizeof(bed5_t *));
    }
    else if ((*lp)->nelems >= (*lp)->capacity) {
        if (bs_globals.verbose) fprintf(stderr, "resizing lookup table [capacity %d -> %d]...\n", (int)(*lp)->capacity, (int)(*lp)->capacity * 2);
        (*lp)->capacity *= 2;
        bed5_t** new_elems = malloc(sizeof(bed5_t *) * (*lp)->capacity);
        for (uint32_t idx = 0; idx < (*lp)->nelems; idx++) {
            if (bs_globals.verbose) fprintf(stderr, "copying old element to new at index [%u]...\n", idx + 1);
            bs_initialize_bed5_element((*lp)->elems[idx]->chr,
                                       (*lp)->elems[idx]->start,
                                       (*lp)->elems[idx]->stop,
                                       (*lp)->elems[idx]->id,
                                       &new_elems[idx]);
            if (bs_globals.verbose) fprintf(stderr, "deleting old element at index [%u]...\n", idx + 1);
            bs_delete_bed5_element(&((*lp)->elems[idx]));
        }   
        (*lp)->elems = new_elems;     
    }
    if (bs_globals.verbose) fprintf(stderr, "number of elements now [%d]\n", (*lp)->nelems + 1);
    uint32_t n = (*lp)->nelems;
    (*lp)->elems[n] = e;
    (*lp)->nelems++;
}

void
bs_delete_lookup_avx(lookup_avx_t** l)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting lookup_avx table...\n");

    for (uint32_t idx = 0; idx < (*l)->nelems; idx++) {
        if (bs_globals.verbose) fprintf(stderr, "deleting element at index [%d]...\n", idx + 1);
        bs_delete_bed5_avx_element(&(*l)->elems[idx]);
    }
    free((*l)->elems);
    (*l)->elems = NULL;
    free(*l);
    *l = NULL;

    if (bs_globals.verbose) fprintf(stderr, "deleted lookup_avx table...\n");
}

void
bs_delete_lookup(lookup_t** l)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting lookup table...\n");

    for (uint32_t idx = 0; idx < (*l)->nelems; idx++) {
        if (bs_globals.verbose) fprintf(stderr, "deleting element at index [%d]...\n", idx + 1);
        bs_delete_bed5_element(&(*l)->elems[idx]);
    }
    free((*l)->elems);
    (*l)->elems = NULL;
    free(*l);
    *l = NULL;

    if (bs_globals.verbose) fprintf(stderr, "deleted lookup table...\n");
}

void
bs_delete_bed5_avx_element(bed5_avx_t** e)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting bed5_avx element...\n");

    free((*e)->chr);
    (*e)->chr = NULL;
    free((*e)->id);
    (*e)->id = NULL;
    if ((*e)->signal) {
        bs_delete_signal_avx(&((*e)->signal));
        (*e)->signal = NULL;
    }
    free(*e), *e = NULL;
}

void
bs_delete_bed5_element(bed5_t** e)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting bed5 element...\n");

    free((*e)->chr);
    (*e)->chr = NULL;
    free((*e)->id);
    (*e)->id = NULL;
    if ((*e)->signal) {
        bs_delete_signal(&((*e)->signal));
        (*e)->signal = NULL;
    }
    free(*e), *e = NULL;
}

void
bs_delete_signal_avx(signal_avx_t** s)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting signal_avx...\n");

    free((*s)->data);
    (*s)->data = NULL;
    free(*s);
    *s = NULL;
}

void
bs_delete_signal(signal_t** s)
{
    if (bs_globals.verbose) fprintf(stderr, "deleting signal...\n");

    free((*s)->data);
    (*s)->data = NULL;
    free(*s);
    *s = NULL;
}

void
bs_initialize_globals()
{
    bs_globals.verbose = false;
    bs_globals.input_fn = NULL;
    bs_globals.output_fn = NULL;
    bs_globals.output_stream = NULL;
    bs_globals.avx_enabled = false;
    bs_globals.lookup = NULL;
    bs_globals.lookup_avx = NULL;
}

void
bs_delete_globals()
{
    if (bs_globals.input_fn) {
        free(bs_globals.input_fn);
        bs_globals.input_fn = NULL;
    }
    if (bs_globals.output_stream) {
        fclose(bs_globals.output_stream);
        bs_globals.output_stream = NULL;
    }
    if (bs_globals.output_fn) {
        free(bs_globals.output_fn);
        bs_globals.output_fn = NULL;
    }
    if (bs_globals.lookup) {
        bs_delete_lookup(&bs_globals.lookup);
    }
    if (bs_globals.lookup_avx) {
        bs_delete_lookup_avx(&bs_globals.lookup_avx);
    }
}

void 
bs_initialize_command_line_options(int argc, char** argv)
{
    int bs_client_long_index;
    int bs_client_opt = getopt_long(argc,
                                    argv,
                                    bs_client_opt_string,
                                    bs_client_long_options,
                                    &bs_client_long_index);

    opterr = 0;

    while (bs_client_opt != -1) {
        switch (bs_client_opt) {
        case 'a':
            bs_globals.avx_enabled = true;
            break;
        case 'i':
            if (!optarg) {
                exit(EXIT_FAILURE);
            }
            bs_globals.input_fn = malloc(strlen(optarg) + 1);
            memcpy(bs_globals.input_fn, optarg, strlen(optarg) + 1);
            break;
        case 'o':
            if (!optarg) {
                exit(EXIT_FAILURE);
            }
            bs_globals.output_fn = malloc(strlen(optarg) + 1);
            memcpy(bs_globals.output_fn, optarg, strlen(optarg) + 1);
            break;
        case 'v':
            bs_globals.verbose = true;
            break;
        case 'h':
        case '?':
            bs_print_usage(stdout);
            exit(EXIT_SUCCESS);
        default:
            break;
        }        
        bs_client_opt = getopt_long(argc,
                                    argv,
                                    bs_client_opt_string,
                                    bs_client_long_options,
                                    &bs_client_long_index);
    }

    if (!bs_globals.output_fn) {
        bs_globals.output_stream = stdout;
    }

    if (!bs_globals.input_fn) {
        bs_print_usage(stdout);
        exit(EXIT_FAILURE);
    }
}

void 
bs_print_usage(FILE* os) 
{
    fprintf(os, "Usage:\n\t$ pearson-test --input <input.bed5> [--output <output.bs>]\n");
}

score_t
bs_pearson_r_via_signal_avx_t(signal_avx_t* a, signal_avx_t* b)
{
    if (a->n != b->n) {
        fprintf(stderr, "Error: Vectors being correlated are of unequal length!\n");
        exit(EXIT_FAILURE);
    }
    if ((a->sd == 0.0f) || (b->sd == 0.0f)) {
        return NAN;
    }
    score_t s = 0.0f;
    
    score_t sum8[8] = {0};
    __m256 a_vec_mean = _mm256_set1_ps(a->mean);
    __m256 b_vec_mean = _mm256_set1_ps(b->mean);
    __m256 sum = _mm256_set1_ps(0.0);
    __m256 a_vec;
    __m256 b_vec;
    score_t* a_ptr = a->data;
    score_t* b_ptr = b->data;
    uint32_t idx = 0;

    for (; idx < (a->n - 8); idx += 8) {
        a_vec = _mm256_load_ps(a_ptr);
        b_vec = _mm256_load_ps(b_ptr);
        sum = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a_vec, a_vec_mean), _mm256_sub_ps(b_vec, b_vec_mean)), sum);
        a_ptr += 8;
        b_ptr += 8;
    }
    uint32_t remaining_bytes = a->n - idx;
    __m256i masks[9];
    masks[0] = _mm256_setr_epi32( 1,  1,  1,  1,  1,  1,  1,  1);
    masks[1] = _mm256_setr_epi32(-1,  1,  1,  1,  1,  1,  1,  1);
    masks[2] = _mm256_setr_epi32(-1, -1,  1,  1,  1,  1,  1,  1);
    masks[3] = _mm256_setr_epi32(-1, -1, -1,  1,  1,  1,  1,  1);
    masks[4] = _mm256_setr_epi32(-1, -1, -1, -1,  1,  1,  1,  1);
    masks[5] = _mm256_setr_epi32(-1, -1, -1, -1, -1,  1,  1,  1);
    masks[6] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1,  1,  1);
    masks[7] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1,  1);
    masks[8] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, -1);

    if (remaining_bytes > 8) {
        fprintf(stderr, "remaining_bytes incorrect [%u]\n", remaining_bytes);
        exit(EXIT_FAILURE);
    }
    __m256i mask = masks[remaining_bytes];
    a_vec = _mm256_maskload_ps(a_ptr, mask);
    b_vec = _mm256_maskload_ps(b_ptr, mask);
    sum = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a_vec, a_vec_mean), _mm256_sub_ps(b_vec, b_vec_mean)), sum);
    sum = _mm256_hadd_ps(sum, sum);
    __m256 sumshuffle = _mm256_permute2f128_ps(sum, sum, 0x1);
    sum = _mm256_add_ps(sum, sumshuffle);
    _mm256_store_ps(sum8, sum);
    s = sum8[0] + sum8[1];

    return s / ((a->n - 1.0f) * a->sd * b->sd);
}

score_t
bs_pearson_r_via_signal_t(signal_t* a, signal_t* b)
{
    if (a->n != b->n) {
        fprintf(stderr, "Error: Vectors being correlated are of unequal length!\n");
        exit(EXIT_FAILURE);
    }
    if ((a->sd == 0.0f) || (b->sd == 0.0f)) {
        return NAN;
    }
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < a->n; idx++)
        s += (a->data[idx] - a->mean) * (b->data[idx] - b->mean);
    return s / ((a->n - 1.0f) * a->sd * b->sd);
}

static inline score_t
bs_mean_signal(score_t* d, uint32_t len)
{
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < len; idx++) {
        s += d[idx];
    }
    return s / len;
}

score_t
bs_sample_sd_signal(score_t* d, uint32_t len, score_t m)
{
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < len; idx++) {
        s += (d[idx] - m) * (d[idx] - m);
    }
    return sqrt(s / (len - 1));
}

inline score_t
bs_truncate_score_to_precision(score_t d, int prec)
{
    score_t factor = (score_t) powf(10, prec);
    return (d < 0) ? (score_t) ceil(d * factor)/factor : (score_t) floor((d + kEpsilon) * factor)/factor;
}

inline byte_t
bs_encode_score_to_byte(score_t d)
{
    if (isnan(d)) {
        return kNANEncodedByte;
    }
    d += (d < 0) ? -kEpsilon : kEpsilon; // jitter is used to deal with interval edges
    d = bs_truncate_score_to_precision(d, 2);
    int encode_d = (int) ((d < 0) ? (ceil(d * 1000.0f)/10.0f + 100) : (floor(d * 1000.0f)/10.0f + 100)) + bs_signbit(-d);
    return (byte_t) encode_d;
}

inline bool
bs_signbit(score_t d)
{
    return (signbit(d) > 0) ? true : false;
}

void
bs_calculate_pearson_scores_via_signal_avx(lookup_avx_t* l)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_calculate_pearson_scores_via_signal_avx()\n");

    for (uint32_t row_idx = 0; row_idx < l->nelems; row_idx++) {
        for (uint32_t col_idx = 0; col_idx < l->nelems; col_idx++) {
            score_t score = bs_pearson_r_via_signal_avx_t(l->elems[row_idx]->signal, l->elems[col_idx]->signal);
            byte_t encoded_score = bs_encode_score_to_byte(score);
            if (bs_globals.verbose) fprintf(stderr, "cell [%05u, %05u] score [%f] encoding [%02x]\n", row_idx, col_idx, score, encoded_score);
            fprintf(stdout, "%c", encoded_score);
        }
    }    
}

void
bs_calculate_pearson_scores_via_signal(lookup_t* l)
{
    if (bs_globals.verbose) fprintf(stderr, "bs_calculate_pearson_scores_via_signal()\n");

    for (uint32_t row_idx = 0; row_idx < l->nelems; row_idx++) {
        for (uint32_t col_idx = 0; col_idx < l->nelems; col_idx++) {
            score_t score = bs_pearson_r_via_signal_t(l->elems[row_idx]->signal, l->elems[col_idx]->signal);
            byte_t encoded_score = bs_encode_score_to_byte(score);
            if (bs_globals.verbose) fprintf(stderr, "cell [%05u, %05u] score [%f] encoding [%02x]\n", row_idx, col_idx, score, encoded_score);
            fprintf(stdout, "%c", encoded_score);
        }
    }    
}

static inline void 
bs_calculate_mean_and_sd_signal_avx(score_t* d, uint32_t n, score_t* m, score_t* sd)
{
    score_t sum8[8] = {0};
    score_t s = 0.0f;
    score_t ss = 0.0f;
    score_t* ptrA = d;
    
    __m256 sum = _mm256_set1_ps(0.0);
    __m256 sumsqr = _mm256_set1_ps(0.0);
    __m256 a;

    for (uint32_t i = 0; i < n; i += 8) {
        a = _mm256_load_ps(ptrA);
        sum = _mm256_add_ps(sum, a);
        a = _mm256_mul_ps(a, a);
        sumsqr = _mm256_add_ps(sumsqr, a);
        ptrA += 8;
    }

    /*
    score_t pOut[8] = {0};
    _mm256_store_ps(pOut, sum);
    fprintf(stderr, "_mm256_store_ps()        [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */
    
    sum = _mm256_hadd_ps(sum, sum);

    /*
    _mm256_store_ps(pOut, sum);
    fprintf(stderr, "_mm256_hadd_ps()         [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    __m256 sumshuffle = _mm256_permute2f128_ps(sum, sum, 0x1);

    /*
    _mm256_store_ps(pOut, sumshuffle);
    fprintf(stderr, "_mm256_permute2f128_ps() [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    sum = _mm256_add_ps(sum, sumshuffle);

    /*
    _mm256_store_ps(pOut, sum);
    fprintf(stderr, "_mm256_add_ps()          [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    _mm256_store_ps(sum8, sum);
    s = sum8[0] + sum8[1];

    /*
    _mm256_store_ps(pOut, sumsqr);
    fprintf(stderr, "_mm256_store_ps()        [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    sumsqr = _mm256_hadd_ps(sumsqr, sumsqr);

    /*
    _mm256_store_ps(pOut, sumsqr);
    fprintf(stderr, "_mm256_hadd_ps()         [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    sumshuffle = _mm256_permute2f128_ps(sumsqr, sumsqr, 0x1);

    /*
    _mm256_store_ps(pOut, sumshuffle);
    fprintf(stderr, "_mm256_permute2f128_ps() [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    sumsqr = _mm256_add_ps(sumsqr, sumshuffle);

    /*
    _mm256_store_ps(pOut, sumsqr);
    fprintf(stderr, "_mm256_add_ps()          [%f %f %f %f %f %f %f %f]\n", 
            pOut[0],
            pOut[1],
            pOut[2],
            pOut[3],
            pOut[4],
            pOut[5],
            pOut[6],
            pOut[7]);
    */

    _mm256_store_ps(sum8, sumsqr);
    ss = sum8[0] + sum8[1];

    *m = s / n;    
    *sd = sqrt((ss - (s * s / n))/(n - 1));
}
