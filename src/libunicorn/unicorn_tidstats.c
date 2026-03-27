#define _XOPEN_SOURCE 700
#include <pthread.h>
#include "unicorn_internal.h"
#include "genesisC.h"

/* ------------------------------------------------------------------ *
 * tid -> taxid cache
 * Built once at startup from the BAM header. Eliminates per-alignment
 * string hash lookups into the acc2tax ensemble map.
 * ------------------------------------------------------------------ */
typedef struct {
    uint32_t *tid2tax;   /* taxid at target rank, indexed by BAM tid  */
    uint8_t  *tid_ok;    /* 1 if this tid has a valid taxid, 0 if not */
    int32_t   n;         /* == hdr->n_targets                         */
} tidcache_t;

static tidcache_t *_tidcache_build(unicorn_t *u, utax_t *utax)
{
    sam_hdr_t *hdr = u->hdr;
    int32_t n = hdr->n_targets;
    tidcache_t *c = calloc(1, sizeof(tidcache_t));
    if (!c) return NULL;
    c->tid2tax = calloc(n, sizeof(uint32_t));
    c->tid_ok  = calloc(n, sizeof(uint8_t));
    c->n       = n;
    if (!c->tid2tax || !c->tid_ok) {
        free(c->tid2tax); free(c->tid_ok); free(c); return NULL;
    }
    const char *rank = utax->rank;
    uint32_t nmiss = 0;
    for (int32_t i = 0; i < n; i++) {
        int absent;
        uint32_t taxid = utax_gettaxid(utax, hdr->target_name[i], &absent);
        if (absent) { nmiss++; continue; }
        if (rank) {
            uint8_t ret;
            taxid = utax_getidatrank(utax, taxid, rank, &ret);
            if (ret) { nmiss++; continue; }
        }
        c->tid2tax[i] = taxid;
        c->tid_ok[i]  = 1;
    }
    if (VERBOSE)
        fprintf(stderr,
                "[libunicorn::%s] tid cache: %d refs, %u missing from taxonomy\n",
                __func__, n, nmiss);
    return c;
}

static void _tidcache_free(tidcache_t *c)
{
    if (!c) return;
    free(c->tid2tax);
    free(c->tid_ok);
    free(c);
}

/* ------------------------------------------------------------------ *
 * refwork_t -- one unit of work: all alignments for one BAM reference.
 * Allocated by the I/O thread, processed by a worker, then freed.
 * ------------------------------------------------------------------ */
typedef struct {
    int32_t      tid;
    uint32_t     taxid;
    uint32_t     reflen;
    ueventq_t    events;
    u64set_t    *readset;
    uint32_t     v_rlen[256];
    float        alnani_mean;
    float        alnani_var;
    float        _MANI;
    float        alnnm_mean;
    float        mdust;
    float        vdust;
    uint64_t     nalns;
    lint2int_t  *camex;      /* NULL when duplicity disabled */
} refwork_t;

/* ------------------------------------------------------------------ *
 * refresult_t -- compact result from worker to merge thread.
 * ------------------------------------------------------------------ */
typedef struct {
    int32_t     tid;
    uint32_t    taxid;
    uint32_t    reflen;
    uint64_t    nalns;
    uint64_t    covbases;
    float       meancov;
    float       meanoncov;
    uint32_t    v_rlen[256];
    float       alnani_mean;
    float       alnani_var;
    float       alnnm_mean;
    float       mdust;
    float       vdust;
    u64set_t   *readset;     /* ownership transferred to merge thread */
    lint2int_t *camex;       /* ownership transferred, or NULL        */
} refresult_t;

/* ------------------------------------------------------------------ *
 * Bounded work queue (one producer / N consumers)
 * ------------------------------------------------------------------ */
typedef struct {
    refwork_t      **ring;
    uint32_t         cap;
    uint32_t         head;
    uint32_t         tail;
    pthread_mutex_t  mu;
    pthread_cond_t   not_empty;
    pthread_cond_t   not_full;
    uint8_t          done;
} workq_t;

static workq_t *_workq_init(uint32_t cap)
{
    workq_t *q = calloc(1, sizeof(workq_t));
    if (!q) return NULL;
    q->ring = calloc(cap, sizeof(refwork_t *));
    if (!q->ring) { free(q); return NULL; }
    q->cap = cap;
    pthread_mutex_init(&q->mu, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
    return q;
}

static void _workq_free(workq_t *q)
{
    if (!q) return;
    pthread_mutex_destroy(&q->mu);
    pthread_cond_destroy(&q->not_empty);
    pthread_cond_destroy(&q->not_full);
    free(q->ring);
    free(q);
}

static int _workq_push(workq_t *q, refwork_t *w)
{
    pthread_mutex_lock(&q->mu);
    while (((q->head + 1) % q->cap) == q->tail && !q->done)
        pthread_cond_wait(&q->not_full, &q->mu);
    if (q->done) { pthread_mutex_unlock(&q->mu); return -1; }
    q->ring[q->head] = w;
    q->head = (q->head + 1) % q->cap;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mu);
    return 0;
}

static refwork_t *_workq_pop(workq_t *q)
{
    pthread_mutex_lock(&q->mu);
    while (q->head == q->tail && !q->done)
        pthread_cond_wait(&q->not_empty, &q->mu);
    if (q->head == q->tail) {
        pthread_mutex_unlock(&q->mu);
        return NULL;
    }
    refwork_t *w = q->ring[q->tail];
    q->tail = (q->tail + 1) % q->cap;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mu);
    return w;
}

static void _workq_seal(workq_t *q)
{
    pthread_mutex_lock(&q->mu);
    q->done = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_cond_broadcast(&q->not_full);
    pthread_mutex_unlock(&q->mu);
}

/* ------------------------------------------------------------------ *
 * Results queue (N producers / 1 consumer)
 * ------------------------------------------------------------------ */
typedef struct {
    refresult_t    **ring;
    uint32_t         cap;
    uint32_t         head;
    uint32_t         tail;
    pthread_mutex_t  mu;
    pthread_cond_t   not_empty;
    pthread_cond_t   not_full;
    uint32_t         n_workers;
    uint32_t         n_done;
} resultq_t;

static resultq_t *_resultq_init(uint32_t cap, uint32_t n_workers)
{
    resultq_t *q = calloc(1, sizeof(resultq_t));
    if (!q) return NULL;
    q->ring = calloc(cap, sizeof(refresult_t *));
    if (!q->ring) { free(q); return NULL; }
    q->cap       = cap;
    q->n_workers = n_workers;
    pthread_mutex_init(&q->mu, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
    return q;
}

static void _resultq_free(resultq_t *q)
{
    if (!q) return;
    pthread_mutex_destroy(&q->mu);
    pthread_cond_destroy(&q->not_empty);
    pthread_cond_destroy(&q->not_full);
    free(q->ring);
    free(q);
}

static void _resultq_push(resultq_t *q, refresult_t *r)
{
    pthread_mutex_lock(&q->mu);
    while (((q->head + 1) % q->cap) == q->tail)
        pthread_cond_wait(&q->not_full, &q->mu);
    q->ring[q->head] = r;
    q->head = (q->head + 1) % q->cap;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mu);
}

static void _resultq_worker_done(resultq_t *q)
{
    pthread_mutex_lock(&q->mu);
    q->n_done++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mu);
}

static refresult_t *_resultq_pop(resultq_t *q)
{
    pthread_mutex_lock(&q->mu);
    while (q->head == q->tail && q->n_done < q->n_workers)
        pthread_cond_wait(&q->not_empty, &q->mu);
    if (q->head == q->tail) {
        pthread_mutex_unlock(&q->mu);
        return NULL;
    }
    refresult_t *r = q->ring[q->tail];
    q->tail = (q->tail + 1) % q->cap;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mu);
    return r;
}

/* ------------------------------------------------------------------ *
 * Worker thread
 * ------------------------------------------------------------------ */
typedef struct {
    workq_t   *wq;
    resultq_t *rq;
} worker_ctx_t;

static refresult_t *_process_ref(refwork_t *w)
{
    refresult_t *r = calloc(1, sizeof(refresult_t));
    if (!r) return NULL;
    r->tid         = w->tid;
    r->taxid       = w->taxid;
    r->reflen      = w->reflen;
    r->nalns       = w->nalns;
    r->alnani_mean = w->alnani_mean;
    r->alnani_var  = w->alnani_var;
    r->alnnm_mean  = w->alnnm_mean;
    r->mdust       = w->mdust;
    r->vdust       = w->vdust;
    r->readset     = w->readset;
    r->camex       = w->camex;
    memcpy(r->v_rlen, w->v_rlen, 256 * sizeof(uint32_t));
    /* Sort events and compute coverage -- always destroy events */
    if (w->events.n) {
        unicorn_sorturange(w->events.n, w->events.a);
        _covstats_t cs = {0};
        _refcoverage(w->events, w->reflen, &cs);
        r->covbases  = cs.covbases;
        r->meancov   = cs.meancov;
        r->meanoncov = cs.meanoncov;
    }
    kv_destroy(w->events); /* always -- kv_init was always called */
    return r;
}

static void *_worker_thread(void *arg)
{
    worker_ctx_t *ctx = (worker_ctx_t *)arg;
    refwork_t *w;
    while ((w = _workq_pop(ctx->wq)) != NULL) {
        refresult_t *r = _process_ref(w);
        free(w);
        if (r) _resultq_push(ctx->rq, r);
    }
    _resultq_worker_done(ctx->rq);
    return NULL;
}

/* ------------------------------------------------------------------ *
 * Merge refresult_t into taxmap (merge thread only -- no locking needed)
 * ------------------------------------------------------------------ */
static void _merge_result(refresult_t *r, taxmap_t *taxmap)
{
    int absent;
    khint_t k = taxmap_get(taxmap, r->taxid);
    if (k == kh_end(taxmap)) {
        taxstat_t zero = {0};
        zero.readset   = u64set_init();
        zero.readl_min = 0xffffffffU;
        k = taxmap_put(taxmap, r->taxid, &absent);
        kh_val(taxmap, k) = zero;
    }
    taxstat_t *t = &kh_val(taxmap, k);

    t->nrefs++;
    t->reflen += r->reflen;
    t->covbases  += r->covbases;
    t->covmean   += r->meancov;   /* averaged over nrefs in _taxmapstats */
    t->meanoncov += r->meanoncov;

    uint64_t prev_nalns = t->nalns;
    t->nalns += r->nalns;

    /* combine Welford means with parallel formula */
    if (r->nalns) {
        double delta;
        delta = r->alnani_mean - t->alnani_mean;
        t->alnani_mean += delta * r->nalns / (double)t->nalns;
        t->alnani_var  += r->alnani_var + delta * delta *
                          ((double)prev_nalns * r->nalns / (double)t->nalns);
        delta = r->alnnm_mean - t->alnnm_mean;
        t->alnnm_mean += delta * r->nalns / (double)t->nalns;
        delta = r->mdust - t->mdust;
        t->mdust += delta * r->nalns / (double)t->nalns;
    }

    for (int i = 0; i < 256; i++)
        t->v_rlen[i] += r->v_rlen[i];

    /* merge readset */
    if (r->readset) {
        khint_t ki;
        kh_foreach(r->readset, ki) {
            int ab;
            uint64_t qid = kh_key(r->readset, ki);
            u64set_put(t->readset, qid, &ab);
        }
        u64set_destroy(r->readset);
        r->readset = NULL;
    }

    /* update readl_min / readl_max */
    for (int i = 0; i < 256; i++) {
        if (r->v_rlen[i]) {
            if ((uint32_t)i < t->readl_min) t->readl_min = i;
            if ((uint32_t)i > t->readl_max) t->readl_max = i;
        }
    }

    /* merge camex if present */
    if (r->camex) {
        if (!t->camex) {
            t->camex = r->camex;
            r->camex = NULL;
        } else {
            khint_t ki;
            kh_foreach(r->camex, ki) {
                uint64_t kmeridx = kh_key(r->camex, ki);
                uint32_t cnt     = kh_val(r->camex, ki);
                int ab;
                khint_t km = lint2int_put(t->camex, kmeridx, &ab);
                if (ab) kh_val(t->camex, km) = cnt;
                else    kh_val(t->camex, km) += cnt;
            }
            lint2int_destroy(r->camex);
            r->camex = NULL;
        }
    }
}

/* ------------------------------------------------------------------ *
 * I/O thread -- reads BAM, groups by tid, pushes refwork_t
 * ------------------------------------------------------------------ */
typedef struct {
    unicorn_t      *u;
    tidcache_t     *cache;
    workq_t        *wq;
    unicorn_stat_t *stats;
    uint8_t         do_duplicity;
    uint8_t         ksize;
    uint64_t        taln;
    uint64_t        kaln;
    uint64_t        nabsent;
} io_ctx_t;

static void *_io_thread(void *arg)
{
    io_ctx_t       *ctx   = (io_ctx_t *)arg;
    unicorn_t      *u     = ctx->u;
    tidcache_t     *cache = ctx->cache;
    unicorn_stat_t *stats = ctx->stats;
    workq_t        *wq    = ctx->wq;
    uint8_t         ksize = ctx->ksize;

    genesis_encoder_t enc = NULL;
    if (ctx->do_duplicity) enc = genesis_encoderinit(ksize);

    bam1_t    *b       = bam_init1();
    refwork_t *cur     = NULL;
    int32_t    cur_tid = -1;
    uint64_t   taln = 0, kaln = 0, nabsent = 0;

    while (sam_read1(u->_FP, u->hdr, b) >= 0) {
        taln++;
        if (_unmapped(b)) continue;
        int32_t tid = b->core.tid;

        /* flush current refwork when reference changes */
        if (tid != cur_tid) {
            if (cur) {
                if (cur->nalns) {
                    _workq_push(wq, cur);
                } else {
                    kv_destroy(cur->events);
                    u64set_destroy(cur->readset);
                    if (cur->camex) lint2int_destroy(cur->camex);
                    free(cur);
                }
            }
            cur     = NULL;
            cur_tid = tid;
        }

        /* taxonomy filter */
        if (!cache->tid_ok[tid]) { nabsent++; continue; }

        /* alignment-level filters */
        if (_reftooshort(u->hdr, tid, stats->minrefl)) continue;
        if (!_ASCHECK(b, stats->minalnas))             continue;
        int32_t dusts = (int)(0.5 + dust(bam_get_seq(b),
                                         b->core.l_qseq, 64, NULL));
        if (dusts > stats->maxdust) continue;

        /* initialise new refwork on first kept alignment for this tid */
        if (!cur) {
            cur = calloc(1, sizeof(refwork_t));
            cur->tid      = tid;
            cur->taxid    = cache->tid2tax[tid];
            cur->reflen   = u->hdr->target_len[tid];
            kv_init(cur->events);
            cur->readset  = u64set_init();
            if (ctx->do_duplicity) cur->camex = lint2int_init();
        }

        uint32_t qlen = b->core.l_qseq;
        kaln++;
        cur->nalns++;

        /* coverage events */
        _urangeevent s = {b->core.pos,   1};
        _urangeevent e = {bam_endpos(b), 0};
        kv_push(_urangeevent, cur->events, s);
        kv_push(_urangeevent, cur->events, e);

        /* read deduplication */
        khint_t q = kh_hash_str(bam_get_qname(b));
        int absent;
        u64set_put(cur->readset, q, &absent);

        if (absent) {
            cur->v_rlen[qlen < 256 ? qlen : 255]++;

            /* kmer duplicity */
            if (ctx->do_duplicity && cur->camex) {
                char seq[256] = {0};
                for (uint32_t i = 0; i < qlen && i < 255; i++)
                    seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
                uint8_t ret;
                for (uint32_t i = 0;
                     i + ksize <= qlen && i + ksize <= 255; i++) {
                    uint64_t kidx = genesis_getcamexidx(enc,
                                                        seq + i,
                                                        ksize, &ret);
                    if (ret) continue;
                    int ab;
                    khint_t km = lint2int_put(cur->camex, kidx, &ab);
                    if (ab) kh_val(cur->camex, km) = 1;
                    else    kh_val(cur->camex, km)++;
                }
            }
        }

        /* ANI / NM / dust -- Welford online, all alignments */
        uint32_t NM;
        float ani = _ANINM(b, &NM);
        uint64_t naln = cur->nalns;
        float delta;
        delta = ani - cur->alnani_mean;
        cur->alnani_mean += delta / naln;
        cur->_MANI       = delta * (ani - cur->alnani_mean);
        cur->alnani_var  = naln > 1 ? cur->_MANI / (naln - 1) : 0.0f;
        delta = (float)NM - cur->alnnm_mean;
        cur->alnnm_mean += delta / naln;
        delta = (float)dusts - cur->mdust;
        cur->mdust += delta / naln;
        cur->vdust  = naln > 1
                    ? delta * ((float)dusts - cur->mdust) / (naln - 1)
                    : 0.0f;
    }

    /* flush last reference */
    if (cur) {
        if (cur->nalns) {
            _workq_push(wq, cur);
        } else {
            kv_destroy(cur->events);
            u64set_destroy(cur->readset);
            if (cur->camex) lint2int_destroy(cur->camex);
            free(cur);
        }
    }

    bam_destroy1(b);
    if (enc) genesis_encoderfree(enc);
    _workq_seal(wq);

    ctx->taln    = taln;
    ctx->kaln    = kaln;
    ctx->nabsent = nabsent;
    return NULL;
}

/* ------------------------------------------------------------------ *
 * _taxmapstats -- finalise per-taxid stats, parallel via kt_forpool
 * ------------------------------------------------------------------ */
typedef struct {
    taxmap_t       *taxmap;
    unicorn_stat_t *stats;
    khint_t        *keys;
    uint32_t        nkeys;
    int32q_t        rmq;
    pthread_mutex_t rmq_mu;
} finalise_ctx_t;

static void _finalise_worker(void *data, long i, int tid)
{
    (void)tid;
    finalise_ctx_t *ctx = (finalise_ctx_t *)data;
    khint_t k = ctx->keys[i];
    if (!kh_exist(ctx->taxmap, k)) return;

    taxstat_t *t    = &kh_val(ctx->taxmap, k);
    uint32_t taxid  = kh_key(ctx->taxmap, k);
    uint32_t nreads = kh_size(t->readset);

    if (nreads < ctx->stats->minnreads ||
        t->alnani_mean < ctx->stats->minmani) {
        pthread_mutex_lock(&ctx->rmq_mu);
        kv_push(int32_t, ctx->rmq, (int32_t)taxid);
        pthread_mutex_unlock(&ctx->rmq_mu);
        return;
    }

    /* finalise coverage averages */
    if (t->nrefs) {
        t->covmean   /= t->nrefs;
        t->meanoncov /= t->nrefs;
    }

    /* median and mode from count histogram */
    t->readl_median = _udCAMEDIAN(t->v_rlen, 256, nreads);
    t->readl_mode   = _udCAMODE(t->v_rlen, 256);

    /* read length mean and variance from histogram */
    t->readl_mean = 0.0f;
    uint32_t total = 0;
    for (int j = 0; j < 256; j++) {
        t->readl_mean += (float)j * t->v_rlen[j];
        total += t->v_rlen[j];
    }
    if (total) t->readl_mean /= total;
    float var = 0.0f;
    for (int j = 0; j < 256; j++) {
        if (!t->v_rlen[j]) continue;
        float d = j - t->readl_mean;
        var += t->v_rlen[j] * d * d;
    }
    t->readl_var = total > 1 ? var / (total - 1) : 0.0f;

    /* duplicity */
    if (t->camex) {
        uint64_t nkmers = 0;
        khint_t ki;
        kh_foreach(t->camex, ki) nkmers += kh_val(t->camex, ki);
        t->duplicity = nkmers
                     ? (float)kh_size(t->camex) / nkmers
                     : 0.0f;
    }
}

static void _taxmapstats(unicorn_stat_t *stats, int threads)
{
    taxmap_t *taxmap = (taxmap_t *)stats->__map;
    if (!kh_size(taxmap)) return;

    /* collect bucket indices for kt_forpool */
    khint_t *keys = malloc(kh_end(taxmap) * sizeof(khint_t));
    uint32_t nkeys = 0;
    khint_t k;
    kh_foreach(taxmap, k) keys[nkeys++] = k;

    finalise_ctx_t ctx = {0};
    ctx.taxmap = taxmap;
    ctx.stats  = stats;
    ctx.keys   = keys;
    ctx.nkeys  = nkeys;
    kv_init(ctx.rmq);
    pthread_mutex_init(&ctx.rmq_mu, NULL);

    void *pool = kt_forpool_init(threads);
    kt_forpool(pool, _finalise_worker, &ctx, nkeys);
    kt_forpool_destroy(pool);

    /* remove filtered taxids and free their resources */
    for (uint32_t i = 0; i < ctx.rmq.n; i++) {
        k = taxmap_get(taxmap, (uint32_t)ctx.rmq.a[i]);
        if (k != kh_end(taxmap)) {
            taxstat_t t = kh_val(taxmap, k);
            u64set_destroy(t.readset);
            if (t.camex) lint2int_destroy(t.camex);
            taxmap_del(taxmap, k);
        }
    }

    pthread_mutex_destroy(&ctx.rmq_mu);
    kv_destroy(ctx.rmq);
    free(keys);
    stats->_nfrefs = kh_size(taxmap);
}

/* ------------------------------------------------------------------ *
 * Public entry point
 * ------------------------------------------------------------------ */
int unicorn_tidstat_compute(unicorn_t *u,
                            unicorn_stat_t *stats,
                            utax_t *utax)
{
    int ret = -1;
    if (!u || !stats || !utax) return ret;

    /* require coordinate-sorted BAM */
    if (!(u->sorted & COORDSORTED)) {
        fprintf(stderr,
                "[libunicorn::%s] ERROR: tidstats requires "
                "coordinate-sorted BAM.\n"
                "  Sort with: samtools sort -o sorted.bam input.bam\n",
                __func__);
        return -2;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* 1. build tid->taxid cache */
    if (VERBOSE)
        fprintf(stderr,
                "[libunicorn::%s] Building tid->taxid cache\n", __func__);
    tidcache_t *cache = _tidcache_build(u, utax);
    if (!cache) return -3;

    /* 2. set up pipeline */
    int threads   = u->threads > 1 ? u->threads : 1;
    int n_workers = threads > 1 ? threads - 1 : 1;
    uint32_t qdepth = (uint32_t)(n_workers * 2 + 4);
    workq_t   *wq = _workq_init(qdepth + 1);
    resultq_t *rq = _resultq_init(qdepth + 1, n_workers);
    if (!wq || !rq) {
        _tidcache_free(cache);
        _workq_free(wq);
        _resultq_free(rq);
        return -3;
    }

    /* 3. start worker threads */
    pthread_t    *workers = calloc(n_workers, sizeof(pthread_t));
    worker_ctx_t  wctx   = { wq, rq };
    for (int i = 0; i < n_workers; i++)
        pthread_create(&workers[i], NULL, _worker_thread, &wctx);

    /* 4. start I/O thread */
    io_ctx_t ioctx = {
        .u            = u,
        .cache        = cache,
        .wq           = wq,
        .stats        = stats,
        .do_duplicity = stats->ksize > 0,
        .ksize        = stats->ksize,
        .taln         = 0,
        .kaln         = 0,
        .nabsent      = 0,
    };
    pthread_t io_tid;
    pthread_create(&io_tid, NULL, _io_thread, &ioctx);

    /* 5. this thread is the merge thread -- drain result queue */
    taxmap_t    *taxmap = (taxmap_t *)stats->__map;
    refresult_t *r;
    while ((r = _resultq_pop(rq)) != NULL) {
        _merge_result(r, taxmap);
        free(r);
    }

    /* 6. join all threads */
    pthread_join(io_tid, NULL);
    for (int i = 0; i < n_workers; i++)
        pthread_join(workers[i], NULL);
    free(workers);

    clock_gettime(CLOCK_MONOTONIC, &t1);

    stats->_nalns  = ioctx.taln;
    stats->_nreads = 0; /* computed below after _taxmapstats */

    if (VERBOSE) {
        uint64_t ns = (t1.tv_sec  - t0.tv_sec)  * 1000000000ULL
                    + (t1.tv_nsec - t0.tv_nsec);
        fprintf(stderr,
                "[libunicorn::%s] BAM pass: %"PRIu64" alns "
                "(%"PRIu64" kept, %"PRIu64" no-tax) in %.2fs\n",
                __func__,
                ioctx.taln, ioctx.kaln, ioctx.nabsent,
                ns / 1e9);
        fprintf(stderr,
                "[libunicorn::%s] %u taxids observed\n",
                __func__, kh_size(taxmap));
    }

    if (!ioctx.kaln) { ret = 0; goto exit; }

    /* 7. finalise per-taxid stats in parallel */
    if (VERBOSE)
        fprintf(stderr,
                "[libunicorn::%s] Finalising taxid stats\n", __func__);
    _taxmapstats(stats, threads);

    /* 8. count total unique reads across all surviving taxids */
    {
        u64set_t *gset = u64set_init();
        int ab;
        khint_t k;
        kh_foreach(taxmap, k) {
            taxstat_t *t = &kh_val(taxmap, k);
            khint_t ki;
            kh_foreach(t->readset, ki)
                u64set_put(gset, kh_key(t->readset, ki), &ab);
        }
        stats->_nreads = kh_size(gset);
        u64set_destroy(gset);
    }

    stats->fc = 1;
    ret = 0;

exit:
    _tidcache_free(cache);
    _workq_free(wq);
    _resultq_free(rq);
    return ret;
}

/* ------------------------------------------------------------------ *
 * Print
 * ------------------------------------------------------------------ */
void unicorn_taxstat_print(const unicorn_t *u,
                           const unicorn_stat_t *stats,
                           FILE *fp,
                           utax_t *utax)
{
    if (!stats || !fp || !u || !utax) return;
    if (!stats->fc) return;
    uint8_t do_dup = stats->ksize > 0;
    fprintf(fp, do_dup ? TIDSTATSTR : TIDSTATSTR_NODUP);
    khint_t k;
    taxmap_t *taxmap = (taxmap_t *)stats->__map;
    kh_foreach(taxmap, k) {
        taxstat_t taxstat = kh_val(taxmap, k);
        float breath    = taxstat.covbases / (double)taxstat.reflen;
        float expbreath = 1.0f - expf(-breath);
        if (do_dup) {
            fprintf(fp, TIDFMTSTR,
                    kh_key(taxmap, k),
                    utax_getname(utax, kh_key(taxmap, k)),
                    taxstat.nrefs,
                    taxstat.reflen,
                    taxstat.nalns,
                    kh_size(taxstat.readset),
                    taxstat.readl_mean,
                    sqrtf(taxstat.readl_var),
                    taxstat.readl_median,
                    taxstat.readl_mode,
                    taxstat.readl_min,
                    taxstat.readl_max,
                    taxstat.alnnm_mean,
                    taxstat.alnani_mean,
                    sqrtf(taxstat.alnani_var),
                    taxstat.covbases,
                    taxstat.covmean,
                    breath,
                    expbreath,
                    breath / expbreath,
                    taxstat.meanoncov,
                    1000.0f * breath,
                    taxstat.duplicity);
        } else {
            fprintf(fp, TIDFMTSTR_NODUP,
                    kh_key(taxmap, k),
                    utax_getname(utax, kh_key(taxmap, k)),
                    taxstat.nrefs,
                    taxstat.reflen,
                    taxstat.nalns,
                    kh_size(taxstat.readset),
                    taxstat.readl_mean,
                    sqrtf(taxstat.readl_var),
                    taxstat.readl_median,
                    taxstat.readl_mode,
                    taxstat.readl_min,
                    taxstat.readl_max,
                    taxstat.alnnm_mean,
                    taxstat.alnani_mean,
                    sqrtf(taxstat.alnani_var),
                    taxstat.covbases,
                    taxstat.covmean,
                    breath,
                    expbreath,
                    breath / expbreath,
                    taxstat.meanoncov,
                    1000.0f * breath);
        }
    }
}
