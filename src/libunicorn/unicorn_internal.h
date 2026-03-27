/*
MIT License

Copyright (c) 2025 GeoGenetics
...
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <inttypes.h>
#include <time.h>

#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "klib/khashl.h"
#include "klib/ksort.h"
#include "klib/kthread.h"
#include "klib/kavl.h"
#include "klib/kvec.h"
typedef struct {
  float  score;
  uint32_t al; //Alignment length
  uint32_t tid;
} alnscore_t;
#define kv_pushq(v, x) do {                                         \
        if ((v).n == (v).m) {                                       \
            (v).m = (v).m? (v).m<<1 : 2;                            \
            (v).a = realloc((v).a, sizeof(bam1_t*) * (v).m);        \
            for (uint32_t i = (v).n; i < (v).m; i++) {              \
                (v).a[i] = bam_init1();                             \
            }                                                       \
        }                                                           \
        (v).a[(v).n] = bam_copy1( (v).a[(v).n], (x) );              \
        (v).n++;                                                    \
} while (0)
#define kv_lastq(v) (v).a[(v).n-1]
typedef kvec_t(bam1_t *)    bamq_t;
typedef kvec_t(float)       floatq_t;
typedef kvec_t(uint32_t)    uint32q_t;
typedef kvec_t(int32_t)     int32q_t;
typedef kvec_t(char *)      charq_t;
typedef kvec_t(char *)      strq_t;
typedef kvec_t(alnscore_t)  alnscoreq_t;
typedef kvec_t(alnscoreq_t) dataq_t;

extern uint8_t VERBOSE;

/*
********************************
 * Rango object for coverage computation.
 * Encodes start e==1 or end e==0 of range.
*/
typedef struct _urangeevent {
    uint64_t pos:63;
    uint8_t e:1;
} _urangeevent;
static inline uint8_t _eventlt(_urangeevent a, _urangeevent b)
{
    if (a.pos != b.pos)
        return a.pos < b.pos;
    return a.e > b.e;
}
typedef kvec_t(_urangeevent) ueventq_t;
void unicorn_sorturange(uint32_t n, _urangeevent *a);

/*Sort values for bam files*/
#define UNSRTED      0x00
#define QUERYSORTED  0x01
#define QUERYGROUPED 0x02
#define COORDSORTED  0x04

typedef struct values_t {
  uint64_t naln;
  uint64_t nfaln;
  uint64_t nread;
  uint64_t nfread;
  uint32_t nref;
  uint32_t nfref;
} values_t;

typedef struct {
  int  argc;
  char **argv;
  int  threads;
  char *ifile;
  char *outbam;
  hts_tpool *p;
  htsFile   *_FP;
  bam_hdr_t *hdr;
  uint8_t sorted;
  uint8_t dcache;
  bam1_t *daln;
  values_t values;
} unicorn_t;

uint8_t unicorn_isqgrouped(unicorn_t *u);
int32_t unicorn_reassignload(unicorn_t *u, alnscoreq_t *q);
int32_t unicorn_alnfiltload(unicorn_t *u, alnscoreq_t *q);
dataq_t *unicorn_qloadqueue(unicorn_t *u, uint64_t *naln);
uint8_t unicorn_rewind(unicorn_t *u);

#define _unmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define _reftooshort(hdr, tid, minref)\
          ((hdr)->target_len[(tid)] < (minref) ? 1 : 0)

KHASHL_MAP_INIT(static,
                int32int64map_t,
                int32int64map,
                uint32_t,
                uint64_t,
                kh_hash_uint32,
                kh_eq_generic)
KHASHL_MAP_INIT(static, lint2int_t, lint2int,
                uint64_t, uint32_t,
                kh_hash_uint64, kh_eq_generic)

KHASHL_SET_INIT(static,
                u64set_t, u64set,
                uint64_t,
                kh_hash_dummy, kh_eq_generic)

typedef struct refstat_t {
  uint32_t     REFLEN;
  uint64_t     REFNALNS;
  float        REFREADE;
  float        REFREADV;
  uint32_t     REFREADD;
  uint32_t     REFREADO;
  float        _M;
  uint32_t     REFREADMIN;
  uint32_t     REFREADMAX;
  u64set_t     *READSET;
  float        REFALNNM;
  float        REFALNANIE;
  float        REFALNANIV;
  float        REFALNANID;
  float        REFALNANIO;
  float        _MANI;
  uint64_t     REFCOVB;
  float        REFMCOV;
  float        REFMONCOV;
  float        REFVONCOV;
  float        REFENTROPY;
  float        REFGINI;
  float        REFNENTROP;
  float        REFNGINI;
  float        tad80;
  float        mdust;
  float        vdust;
  floatq_t     aANI;
  uint32_t     aRLEN[256];
  ueventq_t    aEVENT;
  int32_t      _ntid;
} refstat_t;

KHASHL_MAP_INIT(static,
                refmap_t, refmap,
                int32_t, refstat_t,
                kh_hash_uint32, kh_eq_generic)

typedef struct taxstat_t {
  uint32_t     nrefs;
  uint64_t     reflen;
  uint64_t     nalns;
  float        readl_mean;
  float        readl_var;
  uint32_t     readl_median;
  uint32_t     readl_mode;
  float        _M;
  uint32_t     readl_min;
  uint32_t     readl_max;
  u64set_t     *readset;
  refmap_t     *refmap;
  float        alnnm_mean;
  float        alnani_mean;
  float        alnani_var;
  float        alnani_median;
  float        alnani_mode;
  float        _MANI;
  float        mdust;
  float        vdust;
  float        duplicity;
  lint2int_t   *camex;
  uint64_t     covbases;
  float        covmean;
  float        meanoncov;
  float        varoncov;
  float        coventropy;
  float        covgini;
  float        covnentropy;
  float        covngini;
  float        tad80;
  floatq_t     a_ani;
  uint32_t     v_rlen[256];
  ueventq_t    aEVENT;
  int32_t      _ntid;
} taxstat_t;

KHASHL_MAP_INIT(static,
                taxmap_t, taxmap,
                int32_t, taxstat_t,
                kh_hash_uint32, kh_eq_generic)

KHASHL_MAP_INIT(static,
                floatmap_t, floatmap,
                uint32_t, uint64_t,
                kh_hash_uint32, kh_eq_generic)

typedef struct unicorn_stats_t {
  uint8_t fc: 1;
  void *__map;
  uint8_t mapflg;
  uint64_t _nalns;
  uint64_t _nreads;
  uint64_t _nfreads;
  uint64_t _nfalns;
  uint32_t _nrefs;
  uint32_t _nfrefs;
  float    _mrlen;
  float    _vrlen;
  uint32_t _mdrlen;
  uint32_t _morlen;
  uint32_t _readlc[256];
  floatmap_t *_anihist;
  float   _meanani;
  float   _meannm;
  uint64_t _tlen;
  uint64_t _clen;
  uint32_t minnreads;
  uint32_t minrefl;
  float    minmani;
  int32_t  minalnas;
  int32_t  maxdust;
  uint8_t  ksize;
} unicorn_stat_t;

typedef struct _covstats_t {
  uint64_t covbases;
  float    meancov;
  float    meanoncov;
  float    varoncov;
  float    entropy;
  float    gini;
  float    nentropy;
  float    ngini;
  float    tad80;
} _covstats_t;

#define STATSTR "Id\t"\
                "Length\t"\
                "n_alns\t"\
                "n_reads\t"\
                "m_readl\t"\
                "std_readl\t"\
                "md_readl\t"\
                "mo_readl\t"\
                "readl_min\t"\
                "readl_max\t"\
                "m_alnnm\t"\
                "m_alnani\t"\
                "std_alnani\t"\
                "md_alnani\t"\
                "n_covbases\t"\
                "m_cov\t"\
                "breath_cov\t"\
                "exp_breath\t"\
                "breath_ratio\t"\
                "m_covcovered\t"\
                "std_covcovered\t"\
                "evenness_cov\t"\
                "site_density\t"\
                "entropy\t"\
                "gini\t"\
                "n_entropy\t"\
                "n_gini\t"\
                "tad80\t"\
                "mdust\t"\
                "std_dust\n"

#define STATSTR2 "Id\t"\
                 "taxID\t"\
                 "Length\t"\
                 "n_alns\t"\
                 "n_reads\t"\
                 "m_readl\t"\
                 "std_readl\t"\
                 "md_readl\t"\
                 "mo_readl\t"\
                 "readl_min\t"\
                 "readl_max\t"\
                 "m_alnnm\t"\
                 "m_alnani\t"\
                 "std_alnani\t"\
                 "md_alnani\t"\
                 "n_covbases\t"\
                 "m_cov\t"\
                 "breath_cov\t"\
                 "exp_breath\t"\
                 "breath_ratio\t"\
                 "m_covcovered\t"\
                 "std_covcovered\t"\
                 "evenness_cov\t"\
                 "site_density\t"\
                 "entropy\t"\
                 "gini\t"\
                 "n_entropy\t"\
                 "n_gini\t"\
                 "tad80\n"

/* With duplicity column (used when --duplicity is set) */
#define TIDSTATSTR "#taxid\t"\
                   "name\t"\
                   "num_accessions\t"\
                   "total_length\t"\
                   "num_alns\t"\
                   "num_reads\t"\
                   "mean_readl\t"\
                   "stdev_readl\t"\
                   "median_readl\t"\
                   "mode_readl\t"\
                   "readl_min\t"\
                   "readl_max\t"\
                   "mean_alnnm\t"\
                   "mean_alnani\t"\
                   "stdev_alnani\t"\
                   "num_covbases\t"\
                   "mean_cov\t"\
                   "breath_cov\t"\
                   "exp_breath\t"\
                   "breath_ratio\t"\
                   "mean_covcovered\t"\
                   "site_density\t"\
                   "duplicity\n"

/* Without duplicity column (default) */
#define TIDSTATSTR_NODUP "#taxid\t"\
                         "name\t"\
                         "num_accessions\t"\
                         "total_length\t"\
                         "num_alns\t"\
                         "num_reads\t"\
                         "mean_readl\t"\
                         "stdev_readl\t"\
                         "median_readl\t"\
                         "mode_readl\t"\
                         "readl_min\t"\
                         "readl_max\t"\
                         "mean_alnnm\t"\
                         "mean_alnani\t"\
                         "stdev_alnani\t"\
                         "num_covbases\t"\
                         "mean_cov\t"\
                         "breath_cov\t"\
                         "exp_breath\t"\
                         "breath_ratio\t"\
                         "mean_covcovered\t"\
                         "site_density\n"

/* Format string with duplicity */
#define TIDFMTSTR "%u\t"\
                  "%s\t"\
                  "%u\t"\
                  "%"PRIu64"\t"\
                  "%"PRIu64"\t"\
                  "%u\t"\
                  "%f\t"\
                  "%f\t"\
                  "%u\t"\
                  "%u\t"\
                  "%u\t"\
                  "%u\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\t"\
                  "%"PRIu64"\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\t"\
                  "%f\n"

/* Format string without duplicity */
#define TIDFMTSTR_NODUP "%u\t"\
                        "%s\t"\
                        "%u\t"\
                        "%"PRIu64"\t"\
                        "%"PRIu64"\t"\
                        "%u\t"\
                        "%f\t"\
                        "%f\t"\
                        "%u\t"\
                        "%u\t"\
                        "%u\t"\
                        "%u\t"\
                        "%f\t"\
                        "%f\t"\
                        "%f\t"\
                        "%"PRIu64"\t"\
                        "%f\t"\
                        "%f\t"\
                        "%f\t"\
                        "%f\t"\
                        "%f\t"\
                        "%f\n"

uint8_t _ASCHECK(bam1_t *b, int32_t ms);

uint32_t _udCAMEDIAN(uint32_t *v, uint32_t n, uint32_t mcount);
uint32_t _udCAMODE(uint32_t *v, uint32_t n);
float _ANINM(bam1_t *b, uint32_t *NM);
float _tad80(int32int64map_t *hist);
double _getentropy(const int32int64map_t *hist, uint64_t t, float *_ne);
double _getgini(int32int64map_t *hist, uint64_t t, float m, float *_ng);
uint64_t cov_hist(ueventq_t events,
                  int32int64map_t *hist,
                  uint64_t *_tdepthsum,
                  uint64_t *_sumsqdepth,
                  uint32_t *_maxdepth);
uint64_t _refcoverage(ueventq_t events, uint64_t l, _covstats_t *covstats);

KHASHL_MAP_INIT(static, int2int_t, int2int,
                uint32_t, uint32_t,
                kh_hash_uint32, kh_eq_generic)
KHASHL_MAP_INIT(static, int2chr_t, int2chr,
                uint32_t, char *,
                kh_hash_uint32, kh_eq_generic)
KHASHL_MAP_INIT(static, chr2int_t, chr2int,
                const char *, uint32_t,
                kh_hash_str, kh_eq_str)

typedef struct emap_chr2int_t {
    chr2int_t **maps;
    uint8_t   bits;
    uint64_t  size;
    uint8_t   is_ff;
    char      **keys;
} emap_chr2int_t;

void _echr2intdel(emap_chr2int_t *map);
emap_chr2int_t *_echr2intinit(uint8_t bits, uint8_t is_ff);
int _emapwrite(emap_chr2int_t *map, BGZF *fp);
uint8_t _iskhashfp(BGZF *fp);
emap_chr2int_t *_io_loadkhash(BGZF *fp, int *ret);

typedef struct utupple_t {
  uint32_t taxid;
  const char *rank;
  uint8_t rank_val;
} utuple_t;
KHASHL_MAP_INIT(static, uint2tup_t, uint2tup,
                uint32_t, utuple_t,
                kh_hash_uint32, kh_eq_generic)
KHASHL_SET_INIT(static,
                chrset_t, chrset,
                const char *,
                kh_hash_str, kh_eq_str)
KHASHL_MAP_INIT(static, chr2set_t, chr2set,
                const char *, chrset_t *,
                kh_hash_str, kh_eq_str)

typedef struct nodes_t {
  uint2tup_t *map;
  chr2int_t *levelmap;
} nodes_t;

typedef struct utax_t {
  uint32_t numnodes;
  uint64_t numaccs;
  nodes_t  nodes;
  int2chr_t  *namemap;
  emap_chr2int_t *accmap;
  const char *rank;
} utax_t;

uint32_t utax_gettaxid(utax_t *utax, const char *acc, int *absent);
const char *utax_getname(utax_t *utax, uint32_t taxid);
uint32_t utax_getidatrank(utax_t *utax, uint32_t taxid, const char *rank, uint8_t *ret);
double dust(const uint8_t *seq, int32_t l, int32_t window, int32_t *wCount);
uint64_t _getcovbases(ueventq_t events, uint64_t *_depthsum, uint64_t *_sumsqdepth);
