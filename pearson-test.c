#include "pearson-test.h"

int
main(int argc, char** argv) 
{
    pt_initialize_globals();
    pt_initialize_command_line_options(argc, argv);
    if (!pt_globals.avx_enabled) {
        pt_initialize_lookup_via_signal(pt_globals.input_fn, &pt_globals.lookup);
        //pt_calculate_pearson_scores_via_signal(pt_globals.lookup);
    }
    else {
        pt_initialize_lookup_avx_via_signal_avx(pt_globals.input_fn, &pt_globals.lookup_avx);
        //pt_calculate_pearson_scores_via_signal_avx(pt_globals.lookup_avx);
    }
    pt_delete_globals();

    return EXIT_SUCCESS;
}

void
pt_initialize_lookup_avx_via_signal_avx(const char* ifn, lookup_avx_t** lp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_lookup_avx_via_signal_avx()\n");

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
        pt_initialize_bed5_avx_element(chr_str, start_val, stop_val, id_str, &e);
        pt_push_bed5_avx_element_to_lookup_avx(e, &l);
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
pt_initialize_lookup_via_signal(const char* ifn, lookup_t** lp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_lookup_via_signal()\n");

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
        pt_initialize_bed5_element(chr_str, start_val, stop_val, id_str, &e);
        pt_push_bed5_element_to_lookup(e, &l);
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
pt_initialize_bed5_avx_element(char* chr, uint64_t start, uint64_t stop, char* id, bed5_avx_t** ep)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_bed5_avx_element()\n");

    bed5_avx_t *e = NULL;
    e = malloc(sizeof(bed5_avx_t));
    if (!e) {
        fprintf(stderr, "Error: Could not allocate space for bed5_avx element!\n");
        exit(EXIT_FAILURE);
    }
    e->chr = NULL;
    if (strlen(chr) > 0) {
        e->chr = malloc(sizeof(*chr) * strlen(chr) + 1);
        if (!e->chr) {
            fprintf(stderr, "Error: Could not allocate space for bed5_avx element chromosome!\n");
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
            fprintf(stderr,"Error: Could not allocate space for bed5_avx element id!\n");
            exit(EXIT_FAILURE);
        }
        memcpy(e->id, id, strlen(id) + 1);
    }
    (e->id) ? pt_initialize_signal_avx(e->id, &(e->signal)) : NULL;

    if (!e) {
        fprintf(stderr, "Error: Could not parse BED elements into bed5_avx_t struct\n");
        exit(EXIT_FAILURE);
    }
    
    *ep = e;
}

void
pt_initialize_bed5_element(char* chr, uint64_t start, uint64_t stop, char* id, bed5_t** ep)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_bed5_element()\n");

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
    (e->id) ? pt_initialize_signal(e->id, &(e->signal)) : NULL;

    if (!e) {
        fprintf(stderr, "Error: Could not parse BED elements into bed5_t struct\n");
        exit(EXIT_FAILURE);
    }
    
    *ep = e;
}

void
pt_initialize_signal_avx(char* id, signal_avx_t** sp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_signal_avx()\n");

    signal_avx_t* s = NULL;
    s = malloc(sizeof(signal_avx_t));
    if (!s) {
        fprintf(stderr, "Error: Could not allocate space for signal_avx pointer!\n");
        exit(EXIT_FAILURE);
    }
    s->n_raw = 1;
    s->n = 0;
    s->data = NULL;
    s->mean = NAN;
    s->sd = NAN;

    for (uint32_t idx = 0; idx < strlen(id); idx++) {
        if (id[idx] == kSignalDelim) {
            s->n_raw++;
        }
    }

    /* 
       __m256 can store eight float (score_t) values,
       which is stored in a constant called AVX_FLOAT_N
    */

    uint32_t n_padded = s->n_raw + AVX_FLOAT_N - (s->n_raw % AVX_FLOAT_N);
    score_t* scores = NULL;
    scores = malloc(n_padded * sizeof(*scores));
    if (!scores) {
        fprintf(stderr, "Error: Could not allocate space for scores buffer!\n");
        exit(EXIT_FAILURE);
    }
    /* zero the tail end of the scores array -- any remainder won't affect sum */
    for (uint32_t idx = n_padded - AVX_FLOAT_N; idx < n_padded; idx++) {
        scores[idx] = 0.0f;
    }

    /* allocate space for multiple of eight bytes, with room to spare */
    s->n = (uint32_t)((float)(s->n_raw) / AVX_FLOAT_N) + 1;
    s->data = malloc(sizeof(*s->data) * s->n);
    if (!s->data) {
        fprintf(stderr, "Error: Could not allocate space for signal_avx data pointer!\n");
        exit(EXIT_FAILURE);
    }
    char* start = id;
    char* end = id;
    char entry_buf[ENTRY_MAX_LEN];
    uint32_t entry_idx = 0;
    uint32_t data_idx = 0;
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
        sscanf(entry_buf, "%f", &scores[entry_idx++]);
        if (isnan(scores[entry_idx - 1])) {
            data_contains_nan = true;
        }
        start = end + 1;
    } while (!finished_parsing);

    for (uint32_t k = 0; k < s->n_raw; k += AVX_FLOAT_N) {
        s->data[data_idx++] = _mm256_load_ps(scores + k);
    }
   
    if (!data_contains_nan) {
        s->mean = pt_mean_signal_avx(s->data, n_padded);
        /*
        if (s->n >= 2) {
            s->sd = pt_sample_sd_signal_avx(s->data, s->n, s->mean);
        }
        else {
            fprintf(stderr, "Warning: Vector has one value and therefore cannot have a standard deviation!\n");
        }
        */
    }
    else {
        fprintf(stderr, "signal NAN\n");
        exit(EXIT_FAILURE);
    }

    free(scores);
    *sp = s;
}

void
pt_initialize_signal(char* id, signal_t** sp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_initialize_signal()\n");

    signal_t* s = NULL;
    s = malloc(sizeof(signal_t));
    if (!s) {
        fprintf(stderr, "Error: Could not allocate space for signal pointer!\n");
        exit(EXIT_FAILURE);
    }
    s->n = 1;
    s->data = NULL;
    s->mean = NAN;
    s->sd = NAN;

    for (uint32_t idx = 0; idx < strlen(id); idx++) {
        if (id[idx] == kSignalDelim) {
            s->n++;
        }
    }

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
        s->mean = pt_mean_signal(s->data, s->n);
        /*
        if (s->n >= 2) {
            s->sd = pt_sample_sd_signal(s->data, s->n, s->mean);
        }
        else {
            fprintf(stderr, "Warning: Vector has one value and therefore cannot have a standard deviation!\n");
        }
        */
    }

    *sp = s;
}

void
pt_push_bed5_avx_element_to_lookup_avx(bed5_avx_t* e, lookup_avx_t** lp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_push_bed5_avx_element_to_lookup_avx()\n");

    if ((*lp)->capacity == 0) {
        (*lp)->capacity++;
        (*lp)->elems = malloc(sizeof(bed5_avx_t *));
    }
    else if ((*lp)->nelems >= (*lp)->capacity) {
        if (pt_globals.verbose) fprintf(stderr, "resizing lookup table [capacity %d -> %d]...\n", (int)(*lp)->capacity, (int)(*lp)->capacity * 2);
        (*lp)->capacity *= 2;
        bed5_avx_t** new_elems = malloc(sizeof(bed5_avx_t *) * (*lp)->capacity);
        for (uint32_t idx = 0; idx < (*lp)->nelems; idx++) {
            if (pt_globals.verbose) fprintf(stderr, "copying old element to new at index [%u]...\n", idx + 1);
            pt_initialize_bed5_avx_element((*lp)->elems[idx]->chr,
                                           (*lp)->elems[idx]->start,
                                           (*lp)->elems[idx]->stop,
                                           (*lp)->elems[idx]->id,
                                           &new_elems[idx]);
            if (pt_globals.verbose) fprintf(stderr, "deleting old element at index [%u]...\n", idx + 1);
            pt_delete_bed5_avx_element(&((*lp)->elems[idx]));
        }   
        (*lp)->elems = new_elems;     
    }
    if (pt_globals.verbose) fprintf(stderr, "number of elements now [%d]\n", (*lp)->nelems + 1);
    uint32_t n = (*lp)->nelems;
    (*lp)->elems[n] = e;
    (*lp)->nelems++;
}

void
pt_push_bed5_element_to_lookup(bed5_t* e, lookup_t** lp)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_push_bed5_element_to_lookup()\n");

    if ((*lp)->capacity == 0) {
        (*lp)->capacity++;
        (*lp)->elems = malloc(sizeof(bed5_t *));
    }
    else if ((*lp)->nelems >= (*lp)->capacity) {
        if (pt_globals.verbose) fprintf(stderr, "resizing lookup table [capacity %d -> %d]...\n", (int)(*lp)->capacity, (int)(*lp)->capacity * 2);
        (*lp)->capacity *= 2;
        bed5_t** new_elems = malloc(sizeof(bed5_t *) * (*lp)->capacity);
        for (uint32_t idx = 0; idx < (*lp)->nelems; idx++) {
            if (pt_globals.verbose) fprintf(stderr, "copying old element to new at index [%u]...\n", idx + 1);
            pt_initialize_bed5_element((*lp)->elems[idx]->chr,
                                       (*lp)->elems[idx]->start,
                                       (*lp)->elems[idx]->stop,
                                       (*lp)->elems[idx]->id,
                                       &new_elems[idx]);
            if (pt_globals.verbose) fprintf(stderr, "deleting old element at index [%u]...\n", idx + 1);
            pt_delete_bed5_element(&((*lp)->elems[idx]));
        }   
        (*lp)->elems = new_elems;     
    }
    if (pt_globals.verbose) fprintf(stderr, "number of elements now [%d]\n", (*lp)->nelems + 1);
    uint32_t n = (*lp)->nelems;
    (*lp)->elems[n] = e;
    (*lp)->nelems++;
}

void
pt_delete_lookup_avx(lookup_avx_t** l)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting lookup_avx table...\n");

    for (uint32_t idx = 0; idx < (*l)->nelems; idx++) {
        if (pt_globals.verbose) fprintf(stderr, "deleting element at index [%d]...\n", idx + 1);
        pt_delete_bed5_avx_element(&(*l)->elems[idx]);
    }
    free((*l)->elems);
    (*l)->elems = NULL;
    free(*l);
    *l = NULL;

    if (pt_globals.verbose) fprintf(stderr, "deleted lookup_avx table...\n");
}

void
pt_delete_lookup(lookup_t** l)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting lookup table...\n");

    for (uint32_t idx = 0; idx < (*l)->nelems; idx++) {
        if (pt_globals.verbose) fprintf(stderr, "deleting element at index [%d]...\n", idx + 1);
        pt_delete_bed5_element(&(*l)->elems[idx]);
    }
    free((*l)->elems);
    (*l)->elems = NULL;
    free(*l);
    *l = NULL;

    if (pt_globals.verbose) fprintf(stderr, "deleted lookup table...\n");
}

void
pt_delete_bed5_avx_element(bed5_avx_t** e)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting bed5_avx element...\n");

    free((*e)->chr);
    (*e)->chr = NULL;
    free((*e)->id);
    (*e)->id = NULL;
    if ((*e)->signal) {
        pt_delete_signal_avx(&((*e)->signal));
        (*e)->signal = NULL;
    }
    free(*e), *e = NULL;
}

void
pt_delete_bed5_element(bed5_t** e)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting bed5 element...\n");

    free((*e)->chr);
    (*e)->chr = NULL;
    free((*e)->id);
    (*e)->id = NULL;
    if ((*e)->signal) {
        pt_delete_signal(&((*e)->signal));
        (*e)->signal = NULL;
    }
    free(*e), *e = NULL;
}

void
pt_delete_signal_avx(signal_avx_t** s)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting signal_avx...\n");

    free((*s)->data);
    (*s)->data = NULL;
    free(*s);
    *s = NULL;
}

void
pt_delete_signal(signal_t** s)
{
    if (pt_globals.verbose) fprintf(stderr, "deleting signal...\n");

    free((*s)->data);
    (*s)->data = NULL;
    free(*s);
    *s = NULL;
}

void
pt_initialize_globals()
{
    pt_globals.verbose = false;
    pt_globals.input_fn = NULL;
    pt_globals.output_fn = NULL;
    pt_globals.output_stream = NULL;
    pt_globals.avx_enabled = false;
    pt_globals.lookup = NULL;
    pt_globals.lookup_avx = NULL;
}

void
pt_delete_globals()
{
    if (pt_globals.input_fn) {
        free(pt_globals.input_fn);
        pt_globals.input_fn = NULL;
    }
    if (pt_globals.output_stream) {
        fclose(pt_globals.output_stream);
        pt_globals.output_stream = NULL;
    }
    if (pt_globals.output_fn) {
        free(pt_globals.output_fn);
        pt_globals.output_fn = NULL;
    }
    if (pt_globals.lookup) {
        pt_delete_lookup(&pt_globals.lookup);
    }
    if (pt_globals.lookup_avx) {
        pt_delete_lookup_avx(&pt_globals.lookup_avx);
    }
}

void 
pt_initialize_command_line_options(int argc, char** argv)
{
    int pt_client_long_index;
    int pt_client_opt = getopt_long(argc,
                                    argv,
                                    pt_client_opt_string,
                                    pt_client_long_options,
                                    &pt_client_long_index);

    opterr = 0;

    while (pt_client_opt != -1) {
        switch (pt_client_opt) {
        case 'a':
            pt_globals.avx_enabled = true;
            break;
        case 'i':
            if (!optarg) {
                exit(EXIT_FAILURE);
            }
            pt_globals.input_fn = malloc(strlen(optarg) + 1);
            memcpy(pt_globals.input_fn, optarg, strlen(optarg) + 1);
            break;
        case 'o':
            if (!optarg) {
                exit(EXIT_FAILURE);
            }
            pt_globals.output_fn = malloc(strlen(optarg) + 1);
            memcpy(pt_globals.output_fn, optarg, strlen(optarg) + 1);
            break;
        case 'v':
            pt_globals.verbose = true;
            break;
        case 'h':
        case '?':
            pt_print_usage(stdout);
            exit(EXIT_SUCCESS);
        default:
            break;
        }        
        pt_client_opt = getopt_long(argc,
                                    argv,
                                    pt_client_opt_string,
                                    pt_client_long_options,
                                    &pt_client_long_index);
    }

    if (!pt_globals.output_fn) {
        pt_globals.output_stream = stdout;
    }

    if (!pt_globals.input_fn) {
        pt_print_usage(stdout);
        exit(EXIT_FAILURE);
    }
}

void 
pt_print_usage(FILE* os) 
{
    fprintf(os, "Usage:\n\t$ pearson-test --input <input.bed5> [--output <output.bs>]\n");
}

score_t
pt_pearson_r_via_signal_avx_t(signal_avx_t* a, signal_avx_t* b)
{
    if (a->n != b->n) {
        fprintf(stderr, "Error: Vectors being correlated are of unequal length!\n");
        exit(EXIT_FAILURE);
    }
    if ((a->sd == 0.0f) || (b->sd == 0.0f)) {
        return NAN;
    }
    score_t s = 0.0f;
    /*
    for (uint32_t idx = 0; idx < a->n; idx++)
        s += (a->data[idx] - a->mean) * (b->data[idx] - b->mean);
    */
    //__m256 vector_sum = _mm256_setzero_ps();
    return s / ((a->n - 1.0f) * a->sd * b->sd);
}

score_t
pt_pearson_r_via_signal_t(signal_t* a, signal_t* b)
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
pt_mean_signal_avx(__m256* d, uint32_t len)
{
    score_t tmp[AVX_FLOAT_N];
    score_t s;
    const unsigned int eighthPoints = len / 8;

    __m256 accumulator = _mm256_setzero_ps();

    for (uint32_t idx = 0; idx < eighthPoints; idx++) {
        accumulator = _mm256_add_ps(accumulator, d[idx]);
    }

    _mm256_store_ps(tmp, accumulator);

    s = tmp[0];
    s += tmp[1];
    s += tmp[2];
    s += tmp[3];
    s += tmp[4];
    s += tmp[5];
    s += tmp[6];
    s += tmp[7];

    return s / len;
}

score_t
pt_mean_signal(score_t* d, uint32_t len)
{
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < len; idx++) {
        s += d[idx];
    }
    return s / len;
}

score_t
pt_sample_sd_signal_avx(score_t* d, uint32_t len, score_t m)
{
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < len; idx++) {
        s += (d[idx] - m) * (d[idx] - m);
    }
    return sqrt(s / (len - 1));
}

score_t
pt_sample_sd_signal(score_t* d, uint32_t len, score_t m)
{
    score_t s = 0.0f;
    for (uint32_t idx = 0; idx < len; idx++) {
        s += (d[idx] - m) * (d[idx] - m);
    }
    return sqrt(s / (len - 1));
}

inline score_t
pt_truncate_score_to_precision(score_t d, int prec)
{
    score_t factor = (score_t) powf(10, prec);
    return (d < 0) ? (score_t) ceil(d * factor)/factor : (score_t) floor((d + kEpsilon) * factor)/factor;
}

inline byte_t
pt_encode_score_to_byte(score_t d)
{
    if (isnan(d)) {
        return kNANEncodedByte;
    }
    d += (d < 0) ? -kEpsilon : kEpsilon; // jitter is used to deal with interval edges
    d = pt_truncate_score_to_precision(d, 2);
    int encode_d = (int) ((d < 0) ? (ceil(d * 1000.0f)/10.0f + 100) : (floor(d * 1000.0f)/10.0f + 100)) + pt_signbit(-d);
    return (byte_t) encode_d;
}

inline bool
pt_signbit(score_t d)
{
    return (signbit(d) > 0) ? true : false;
}

void
pt_calculate_pearson_scores_via_signal_avx(lookup_avx_t* l)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_calculate_pearson_scores_via_signal_avx()\n");

    for (uint32_t row_idx = 0; row_idx < l->nelems; row_idx++) {
        for (uint32_t col_idx = 0; col_idx < l->nelems; col_idx++) {
            score_t score = pt_pearson_r_via_signal_avx_t(l->elems[row_idx]->signal, l->elems[col_idx]->signal);
            byte_t encoded_score = pt_encode_score_to_byte(score);
            if (pt_globals.verbose) fprintf(stderr, "cell [%05u, %05u] score [%f] encoding [%02x]\n", row_idx, col_idx, score, encoded_score);
            fprintf(stdout, "%c", encoded_score);
        }
    }    
}

void
pt_calculate_pearson_scores_via_signal(lookup_t* l)
{
    if (pt_globals.verbose) fprintf(stderr, "pt_calculate_pearson_scores_via_signal()\n");

    for (uint32_t row_idx = 0; row_idx < l->nelems; row_idx++) {
        for (uint32_t col_idx = 0; col_idx < l->nelems; col_idx++) {
            score_t score = pt_pearson_r_via_signal_t(l->elems[row_idx]->signal, l->elems[col_idx]->signal);
            byte_t encoded_score = pt_encode_score_to_byte(score);
            if (pt_globals.verbose) fprintf(stderr, "cell [%05u, %05u] score [%f] encoding [%02x]\n", row_idx, col_idx, score, encoded_score);
            fprintf(stdout, "%c", encoded_score);
        }
    }    
}
