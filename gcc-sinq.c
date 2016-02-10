/*
 * The goal of this program is to evaluate the precision of sinq(), after
 * seeing some odd results in another piece of code.
 *
 * It uses MPFR as an approximation to an exact sin(x), and to reduce rounding
 * errors in the statistics.
 */

#include <signal.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <quadmath.h>

#define MPFR_PREC ((mpfr_prec_t)512)
#define NOBJS(x) (sizeof(x)/sizeof(*(x)))

/*
 * Various hardcoded representations of values just above PI/2.
 */
static const float gt_pid2 = 0x1.921fb6p+0;
static const uint32_t gt_pid2_int = 0x3fc90fdb;
static const uint32_t gt_pid2_ls23_int = 13176795;



/*
 * These generate a random single-precision float to use as an input for the
 * tests.
 *
 * Generate an integer bitpattern which, when transferred to a 32-bit float,
 * returns a value in +0 <= x < pi/2
 */
static float getrandom_1( gmp_randstate_t state )
{
    union { uint32_t ui; float f; } x;
    x.ui = gmp_urandomm_ui( state, gt_pid2_int );
    return x.f;
}
/*
 * Generate an integer, 0 <= i < 2^23*pi/2, which, when divided by 2^23, turns
 * into a uniformly distributed single-precision float in +0 <= x < pi/2. The
 * division by 2^23 is exact, as is the conversion from integer to float.
 */
static float getrandom_2( gmp_randstate_t state )
{
    float f = gmp_urandomm_ui( state, gt_pid2_ls23_int );
    return f/(1LU<<23);
}
/*
 * Generate random numbers over the whole span of a 32-bit floating point
 * number, excluding NaN and infinites but including denormalized.
 */
static float getrandom_3( gmp_randstate_t state )
{
    union { uint32_t ui; float f; } x;
    do x.ui = gmp_urandomb_ui( state, 32 );
    while( (x.ui & 0x7F800000) == 0x7F800000 );
    return x.f;
}

static const struct {
    const char *name;
    float (* const code)( gmp_randstate_t );
} distribution[] = {
    { "+0 <= x < PI/2, non-uniform", getrandom_1 },
    { "+0 <= x < PI/2, uniform", getrandom_2 },
    { "all floats", getrandom_3 },
};




/*
 * These functions generate a sin(x) from a float input.
 */
static void getsin_flt( mpfr_t y, float phase )
{
    mpfr_set_flt( y, sinf(phase), MPFR_RNDN );
}
static void getsin_d( mpfr_t y, float phase )
{
    mpfr_set_d( y, sin((double)phase), MPFR_RNDN );
}
static void getsin_ld( mpfr_t y, float phase )
{
    mpfr_set_ld( y, sinl((long double)phase), MPFR_RNDN );
}
static void getsin_q( mpfr_t y, float phase )
{
    char s[64];
    quadmath_snprintf( s, sizeof(s), "%Qa", sinq((__float128)phase ) );
    mpfr_set_str( y, s, 0, MPFR_RNDN );
}

static const struct {
    const char *name;
    void (* const code)( mpfr_t, float );
} generator[] = {
    { "float", getsin_flt },
    { "double", getsin_d },
    { "long double", getsin_ld },
    { "__float128", getsin_q },
};






/*
 * Analyzes the error for one value. The inputs are arbitrary precision MPFR
 * types, and a state to track statistics.
 */

typedef struct stats_s {
    const char *distribution, *generator;
    unsigned long n;
    mpfr_t mean, m2;
} stats_t;

static void stats_init( stats_t *stats, const char *dist, const char *gen )
{
    stats->distribution = dist;
    stats->generator = gen;
    stats->n = 0;
    mpfr_inits2( MPFR_PREC, stats->mean, stats->m2, (mpfr_ptr)0 );
    mpfr_set_ui( stats->mean, 0, MPFR_RNDN );
    mpfr_set_ui( stats->m2, 0, MPFR_RNDN );
}

static void stats_clear( stats_t *stats )
{
    mpfr_clears( stats->mean, stats->m2, (mpfr_ptr)0 );
}

static void stats_add( stats_t *stats, const mpfr_t y, const mpfr_t ref )
{
    MPFR_DECL_INIT( t, MPFR_PREC );

    stats->n++;

    /*
     * Calculate the relative change, reldiff = (y-ref)/|ref|
     */
    MPFR_DECL_INIT( diff, MPFR_PREC );
    MPFR_DECL_INIT( reldiff, MPFR_PREC );

    mpfr_sub( diff, y, ref, MPFR_RNDN );
    mpfr_abs( reldiff, ref, MPFR_RNDN );
    mpfr_div( reldiff, diff, reldiff, MPFR_RNDN );

    /*
     * Do some statistics on this.
     */
    MPFR_DECL_INIT( delta, MPFR_PREC );

    mpfr_sub( delta, reldiff, stats->mean, MPFR_RNDN );
    mpfr_div_ui( t, delta, stats->n, MPFR_RNDN );
    mpfr_add( stats->mean, stats->mean, t, MPFR_RNDN );
    mpfr_sub( t, reldiff, stats->mean, MPFR_RNDN );
    mpfr_fma( stats->m2, delta, t, stats->m2, MPFR_RNDN );
}

void stats_print( const stats_t *stats, FILE *fh )
{
    mpfr_fprintf( fh, "#   Distribution: \"%s\"   Generator: \"%s\"\n",
            stats->distribution, stats->generator );
    mpfr_fprintf( fh, "Samples: %lu\n", stats->n );
    mpfr_fprintf( fh, "Relative difference mean: %.10Re\n", stats->mean );

    MPFR_DECL_INIT( variance, MPFR_PREC );
    MPFR_DECL_INIT( stddev, MPFR_PREC );

    if( stats->n > 1 )
    {
        mpfr_div_ui( variance, stats->m2, stats->n-1, MPFR_RNDN );
        mpfr_sqrt( stddev, variance, MPFR_RNDN );
    }

    mpfr_fprintf( fh, "Relative difference variance: %.10Re\n", variance );
    mpfr_fprintf( fh, "Relative difference standard deviation: %.10Re\n",
            stddev );
    mpfr_fprintf( fh, "\n" );
}






int main( void )
{
    /*
     * Abort nicely with CTRL-C and print with HUP.
     */
    volatile bool running = true;
    volatile bool print = false;
    signal( SIGINT, ({ void f(int x){x=x;running=false;};f;}) );
    signal( SIGHUP, ({ void f(int x){x=x;print=true;};f;}) );
    
    /*
     * Collect different statistics for each combination.
     */
    stats_t stats[NOBJS(distribution)][NOBJS(generator)];
    for( size_t dist=0; dist<NOBJS(distribution); dist++ )
        for( size_t gen=0; gen<NOBJS(generator); gen++ )
            stats_init( &stats[dist][gen], distribution[dist].name,
                    generator[gen].name );

    /*
     * Use GMP for random numbers.
     */
    gmp_randstate_t randstate;
    gmp_randinit_mt( randstate );
    gmp_randseed_ui( randstate, 1111 );

    MPFR_DECL_INIT( mpfrtemp, 128+16 );
    MPFR_DECL_INIT( refsin, MPFR_PREC );

    while( running )
    {
        /*
         * Run all combinations of distributions and generators.
         */
        for( size_t dist=0; dist<NOBJS(distribution); dist++ )
        {
            /*
             * Get a random phase and its "reference" sine.
             */
            float phase = distribution[dist].code( randstate );
            mpfr_set_flt( mpfrtemp, phase, MPFR_RNDN );
            mpfr_sin( refsin, mpfrtemp, MPFR_RNDN );

            /*
             * Compare it with the generators.
             */
            for( size_t gen=0; gen<NOBJS(generator); gen++ )
            {
                generator[gen].code( mpfrtemp, phase );
                stats_add( &stats[dist][gen], mpfrtemp, refsin );
            }
        }

        /*
         * Print stats once if HUP is received.
         */
        for( ; print; print=false )
            for( size_t dist=0; dist<NOBJS(distribution); dist++ )
                for( size_t gen=0; gen<NOBJS(generator); gen++ )
                    stats_print( &stats[dist][gen], stdout );
    }

    /*
     * Print the final result.
     */
    for( size_t dist=0; dist<NOBJS(distribution); dist++ )
        for( size_t gen=0; gen<NOBJS(generator); gen++ )
            stats_print( &stats[dist][gen], stdout );

    /*
     * Clean up.
     */
    for( size_t dist=0; dist<NOBJS(distribution); dist++ )
        for( size_t gen=0; gen<NOBJS(generator); gen++ )
            stats_clear( &stats[dist][gen] );
    gmp_randclear( randstate );
    mpfr_free_cache();

    return 0;
}
