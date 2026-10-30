#define _GNU_SOURCE
/*
 * solar_system.c  –  3-D Solar System Simulation
 *
 * Physics  : Full N-body gravity (Sun + 8 planets, all mutual perturbations)
 *            integrated with a 4th-order Yoshida symplectic composition of
 *            kick-drift-kick leapfrog steps (dt = 2 h).  Plain 2nd-order
 *            leapfrog is selectable at run time (key I) for comparison.
 *            Optional 1PN general-relativistic correction (key G) — shows
 *            Mercury's 43″/century perihelion advance live.
 *            Optional swarm of massless test particles: main asteroid belt
 *            and Kuiper belt (key A).
 *
 * Rendering: Software 3-D projection onto an SDL2 2-D renderer at the
 *            display's native resolution (crisp at any size / DPR).
 *            Arcball rotation, pan, zoom, follow-the-body focus mode.
 *
 * Units    : AU (astronomical units), solar masses, Julian years.
 *            GM☉ = k²·365.25² = 39.47693 AU³/yr² (≈4π²);  c = 63241.077 AU/yr.
 *
 * Builds   : one source for both targets — see Makefile.
 *              make            native binary (resizable, HiDPI, video capture)
 *              make wasm       WebAssembly via Emscripten (single file)
 *              make test       physics self-tests (solar_system --test)
 *            Platform differences are confined to #ifdef __EMSCRIPTEN__:
 *            logical window size, font location, main-loop driver, resize
 *            source, and the ffmpeg recorder (native only).
 *
 * SDL3 port: Replace SDL_RenderDrawLineF/LinesF → SDL_RenderLine/Lines and
 *            the init/destroy boilerplate.  Everything else is identical.
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#ifdef __EMSCRIPTEN__
#  include <emscripten/emscripten.h>
#endif

/* Exact J2000 state vectors from JPL Horizons, if tests/fetch_j2000_state.sh
   has been run (the Makefile defines HAVE_HORIZONS_J2000 when the header
   exists).  Otherwise the approximate mean-elements table below is used.  */
#ifdef HAVE_HORIZONS_J2000
#  include "horizons_j2000.h"
#  define IC_SOURCE "JPL Horizons vectors"
#else
#  define IC_SOURCE "JPL mean elements"
#endif

/* ── constants ─────────────────────────────────────────────────────────── */
/* GM☉ in AU³/yr² with the Julian year (365.25 d) as the time unit:
   k²·365.25², k = 0.01720209895 AU^1.5/day (Gaussian constant; equals
   DE441's GM☉ = 1.32712440041e20 m³/s² with the IAU 2012 AU).  The
   textbook G = 4π² is exact only in *Gaussian* years (365.2569 d) and
   would make every mean motion 1.9e-5 too fast against Horizons.       */
#define GM_SUN      39.476926421373    /* AU³/yr²  (= 4π² × 0.9999622)   */
#define G_AU        GM_SUN             /* G with M☉ = 1                  */
#define C_AU_YR     63241.077               /* speed of light, AU/yr        */
#define N_BODIES    9                       /* Sun + 8 planets              */
#define TRAIL_LEN   21000                   /* trail history points         */
#define DT_HOURS    2.0                     /* integration step             */
#define JD_J2000    2451545.0               /* 2000-01-01 12:00 TT          */

/* Logical layout size.  Everything is laid out in this coordinate space and
   mapped to the display's physical pixels by g_S/g_ox/g_oy (see "device
   transform" below), so the output is crisp at any window size or DPR.   */
#ifdef __EMSCRIPTEN__
#  define WIN_W     1280                    /* fits a browser/itch.io embed */
#  define WIN_H     720
#  define CAM_SCALE0 55.0                   /* initial zoom: inner system   */
#else
#  define WIN_W     1600
#  define WIN_H     900
#  define CAM_SCALE0 60.0
#endif
#define PANEL_W     360                     /* right-side info panel width  */
#define SIM_W       (WIN_W - PANEL_W)

/* Directory the Makefile installs the bundled font to (make install).
   Searched before the system font paths so an installed binary is
   self-sufficient.  Overridable: -DSOLAR_DATADIR='"/some/dir"'.           */
#ifndef SOLAR_DATADIR
#  define SOLAR_DATADIR "."
#endif

/* Test particles: integrated as massless bodies under the 9 massive ones,
   with their own (longer) leapfrog step.  Cost is O(N_AST × N_BODIES).    */
#define N_AST_MAIN   1500                   /* main belt, 2.1–3.3 AU        */
#define N_AST_KUIPER  500                   /* classical Kuiper belt, 39–48 */
#define N_AST        (N_AST_MAIN + N_AST_KUIPER)
#define AST_SUB      12                     /* asteroid step = 12 × dt = 1 d */

/* ── body data ─────────────────────────────────────────────────────────── */
typedef struct {
    const char *name;
    double mass;            /* solar masses                                 */
    double a;               /* semi-major axis [AU]                         */
    double ecc;             /* eccentricity                                 */
    double inc;             /* inclination [deg]                            */
    double omega;           /* longitude of ascending node Ω [deg]          */
    double w;               /* argument of perihelion ω [deg]               */
    double M0;              /* mean anomaly at J2000 [deg]                  */
    double radius_km;       /* physical radius [km]  (display only)         */
    double period_yr;       /* sidereal period [yr]  (display only)         */
    double v_orb_km_s;      /* mean orbital speed [km/s] (display only)     */
    Uint8  r, g, b;         /* display colour                               */
    double draw_r;          /* on-screen radius [logical px]                */
} BodyDef;

/* All orbital elements: JPL SSD Table 1, J2000.0 epoch.  Masses in M☉.
   ω = long.peri − long.node,  M0 = mean.lon − long.peri  (both mod 360°) */
static const BodyDef BODIES[N_BODIES] = {
/* Masses: DE441 GM ratios to the Sun — planetary *systems* (Earth = Earth–Moon
   barycentre), matching the barycentric Horizons vectors the model starts from.
    name         mass(M☉)       a(AU)        ecc         inc(°)     Ω(°)        ω(°)       M0(°)      R(km)      T(yr)   v(km/s) R    G    B   dr */
  { "Sun",    1.0,           0.0,        0.0,        0.0,       0.0,        0.0,        0.0,   695700.0,    0.0,    0.0, 255,220, 50, 14.0 },
  { "Mercury",1.66012083e-7, 0.38709927, 0.20563593, 7.00498,  48.33077,   29.12703,  174.79,   2439.7,   0.2408,  47.4, 180,180,180,  3.0 },
  { "Venus",  2.44783828e-6, 0.72333566, 0.00677672, 3.39468,  76.67984,   54.92262,   50.38,   6051.8,   0.6152,  35.0, 230,180, 80,  5.0 },
  { "Earth",  3.04043264e-6, 1.00000261, 0.01671123, 0.00000,   0.0,       102.93768,  357.53,  6371.0,   1.0000,  29.8,  80,130,255,  5.0 },
  { "Mars",   3.22715608e-7, 1.52371034, 0.09339410, 1.84969,  49.55954,  286.49683,   19.39,   3389.5,   1.8809,  24.1, 200, 80, 50,  4.0 },
  { "Jupiter",9.54791910e-4, 5.20288700, 0.04838624, 1.30440, 100.47391,  274.25457,   19.65,  69911.0,  11.8618,  13.1, 220,160, 90, 11.0 },
  { "Saturn", 2.85885670e-4, 9.53667594, 0.05386179, 2.48599, 113.66242,  338.93646,  317.36,  58232.0,  29.4571,   9.7, 210,190,130, 10.0 },
  { "Uranus", 4.36624961e-5,19.18916464, 0.04725744, 0.77264,  74.01693,   96.93735,  142.28,  25362.0,  84.0107,   6.8, 150,220,220,  7.0 },
  { "Neptune",5.15138377e-5,30.06992276, 0.00859048, 1.77004, 131.78423,  273.18054,  259.92,  24622.0, 164.7945,   5.4,  80,100,255,  7.0 },
};

/* ── vector helpers ─────────────────────────────────────────────────────── */
typedef struct { double x, y, z; } Vec3;

static inline Vec3 v3add(Vec3 a, Vec3 b)    { return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 v3sub(Vec3 a, Vec3 b)    { return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 v3scale(Vec3 a, double s){ return (Vec3){a.x*s,   a.y*s,   a.z*s  }; }
static inline double v3dot(Vec3 a, Vec3 b)  { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double v3len(Vec3 a)          { return sqrt(v3dot(a,a)); }
static inline Vec3 v3norm(Vec3 a)           { double l=v3len(a); return l>0?v3scale(a,1.0/l):(Vec3){0,0,0}; }
static inline Vec3 v3cross(Vec3 a, Vec3 b){
    return (Vec3){a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

/* ── rotation matrix (row-major 3×3) ────────────────────────────────────── */
typedef struct { double m[3][3]; } Mat3;

static Mat3 mat3_identity(void){
    Mat3 r; memset(&r,0,sizeof(r));
    r.m[0][0]=r.m[1][1]=r.m[2][2]=1.0;
    return r;
}
static Vec3 mat3_mul_vec(const Mat3 *M, Vec3 v){
    return (Vec3){
        M->m[0][0]*v.x + M->m[0][1]*v.y + M->m[0][2]*v.z,
        M->m[1][0]*v.x + M->m[1][1]*v.y + M->m[1][2]*v.z,
        M->m[2][0]*v.x + M->m[2][1]*v.y + M->m[2][2]*v.z
    };
}
static Mat3 mat3_mul(const Mat3 *A, const Mat3 *B){
    Mat3 C; memset(&C,0,sizeof(C));
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++)
        C.m[i][j] += A->m[i][k]*B->m[k][j];
    return C;
}
/* rotation about arbitrary axis (Rodrigues) */
static Mat3 mat3_rot_axis(Vec3 axis, double angle){
    Vec3 u = v3norm(axis);
    double c=cos(angle), s=sin(angle), t=1-c;
    double x=u.x,y=u.y,z=u.z;
    return (Mat3){{
        {t*x*x+c,   t*x*y-s*z, t*x*z+s*y},
        {t*x*y+s*z, t*y*y+c,   t*y*z-s*x},
        {t*x*z-s*y, t*y*z+s*x, t*z*z+c  }
    }};
}

/* ======================================================================
   Physics
   ====================================================================== */
enum { INTEG_LEAPFROG = 0, INTEG_YOSHIDA4 = 1 };

typedef struct {
    Vec3   pos[N_BODIES];
    Vec3   vel[N_BODIES];     /* synchronised with pos (velocity Verlet)    */
    Vec3   acc[N_BODIES];     /* cached acceleration at pos                 */
    double mass[N_BODIES];

    /* trail ring buffers */
    Vec3  trail[N_BODIES][TRAIL_LEN];
    int   trail_head[N_BODIES];
    int   trail_count[N_BODIES];
    int   trail_cap[N_BODIES];  /* per-body capacity (≤ TRAIL_LEN)        */

    double time_yr;           /* elapsed simulation time since J2000 [yr]   */
    double dt;                /* time step [yr]; negative = time reversed   */
    int    integrator;        /* INTEG_*                                    */
    int    gr;                /* 1PN Sun–planet correction on/off           */

    double E0, L0;            /* initial energy and |angular momentum|      */
} Sim;

/* ── solve Kepler's equation via Newton iteration ───────────────────────── */
static double solve_kepler(double M, double e){
    double E = M;
    for(int i=0;i<50;i++){
        double dE = (M - E + e*sin(E))/(1.0 - e*cos(E));
        E += dE;
        if(fabs(dE)<1e-12) break;
    }
    return E;
}

/* Orbital elements → Cartesian state in the ecliptic J2000 frame.
   Angles in radians; mu = G·(M_central + m).  The in-plane velocity uses the
   exact Keplerian expression v = n·a²/r with n = √(mu/a³) — consistent with
   the dynamics actually integrated, so vis-viva holds to round-off at t=0. */
static void kepler_to_cartesian(double a, double e, double inc, double Om,
                                double w, double M0, double mu,
                                Vec3 *pos, Vec3 *vel)
{
    double E   = solve_kepler(M0, e);
    double nu  = 2.0*atan2(sqrt(1+e)*sin(E/2.0), sqrt(1-e)*cos(E/2.0));
    double r   = a*(1 - e*cos(E));
    double px  = r*cos(nu), py = r*sin(nu);

    double n   = sqrt(mu/(a*a*a));                   /* mean motion [rad/yr] */
    /* Correct Keplerian formula (Bate, Mueller & White): vfac = n·a²/r.
       NOT a·n/√(1−e²), which only holds at the semi-latus rectum.         */
    double vfac = n*a*a / r;
    double vx  = -vfac*sin(E);
    double vy  =  vfac*sqrt(1-e*e)*cos(E);

    /* perifocal → ecliptic: R3(−Ω)·R1(−i)·R3(−ω) */
    double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);
    double Qxx = cO*cw - sO*sw*ci,  Qxy = -cO*sw - sO*cw*ci;
    double Qyx = sO*cw + cO*sw*ci,  Qyy = -sO*sw + cO*cw*ci;
    double Qzx = sw*si,              Qzy =  cw*si;

    *pos = (Vec3){ Qxx*px + Qxy*py, Qyx*px + Qyy*py, Qzx*px + Qzy*py };
    *vel = (Vec3){ Qxx*vx + Qxy*vy, Qyx*vx + Qyy*vy, Qzx*vx + Qzy*vy };
}

#ifndef HAVE_HORIZONS_J2000
static void elements_to_cartesian(const BodyDef *bd, Vec3 *pos, Vec3 *vel){
    if(bd->a == 0.0){ *pos=(Vec3){0,0,0}; *vel=(Vec3){0,0,0}; return; }
    double d = M_PI/180.0;
    kepler_to_cartesian(bd->a, bd->ecc, bd->inc*d, bd->omega*d, bd->w*d,
                        bd->M0*d, G_AU*(BODIES[0].mass + bd->mass), pos, vel);
}
#endif

/* ── accelerations: Newtonian pairs + optional 1PN Sun–planet term ─────── */
static void compute_accel(const Vec3 *pos, const Vec3 *vel, const double *mass,
                          Vec3 *acc, int gr)
{
    for(int i=0;i<N_BODIES;i++) acc[i]=(Vec3){0,0,0};
    for(int i=0;i<N_BODIES;i++){
        for(int j=i+1;j<N_BODIES;j++){
            Vec3 r  = v3sub(pos[j], pos[i]);
            double d2 = v3dot(r,r) + 1e-10;   /* softening                */
            double d3 = d2*sqrt(d2);
            double fi = G_AU*mass[j]/d3;
            double fj = G_AU*mass[i]/d3;
            acc[i] = v3add(acc[i], v3scale(r,  fi));
            acc[j] = v3add(acc[j], v3scale(r, -fj));
        }
    }
    if(!gr) return;
    /* First post-Newtonian correction for each planet in the Sun's field
       (Einstein–Infeld–Hoffmann, two-body dominant term):
         a_1PN = μ/(c²r³) · [ (4μ/r − v²) r + 4 (r·v) v ]
       This is what produces Mercury's extra 42.98″/century.  Planet–planet
       1PN terms are ~10⁻⁴ of this and omitted.                           */
    for(int i=1;i<N_BODIES;i++){
        Vec3 r  = v3sub(pos[i], pos[0]);
        Vec3 v  = v3sub(vel[i], vel[0]);
        double mu = G_AU*(mass[0]+mass[i]);
        double r2 = v3dot(r,r), rr = sqrt(r2);
        double v2 = v3dot(v,v), rv = v3dot(r,v);
        double k  = mu/(C_AU_YR*C_AU_YR*r2*rr);
        Vec3 ag   = v3scale(v3add(v3scale(r, 4.0*mu/rr - v2), v3scale(v, 4.0*rv)), k);
        acc[i] = v3add(acc[i], ag);
        acc[0] = v3sub(acc[0], v3scale(ag, mass[i]/mass[0]));   /* reaction */
    }
}

/* ── one kick-drift-kick step of length h (velocity Verlet) ─────────────── */
static void kdk_step(Sim *sim, double h){
    double hh = 0.5*h;
    for(int i=0;i<N_BODIES;i++)                       /* kick:  v += a·h/2 */
        sim->vel[i] = v3add(sim->vel[i], v3scale(sim->acc[i], hh));
    for(int i=0;i<N_BODIES;i++)                       /* drift: x += v·h   */
        sim->pos[i] = v3add(sim->pos[i], v3scale(sim->vel[i], h));
    compute_accel(sim->pos, sim->vel, sim->mass, sim->acc, sim->gr);
    for(int i=0;i<N_BODIES;i++)                       /* kick:  v += a·h/2 */
        sim->vel[i] = v3add(sim->vel[i], v3scale(sim->acc[i], hh));
}

/* ── advance by one time step sim->dt ───────────────────────────────────── */
static void sim_step(Sim *sim){
    if(sim->integrator == INTEG_YOSHIDA4){
        /* Yoshida (1990) 4th-order composition: three leapfrog sub-steps of
           length w1·dt, w0·dt, w1·dt with w1 = 1/(2−2^{1/3}), w0 = 1−2·w1.
           Error ∝ dt⁴: at dt = 2 h Mercury's numerical perihelion drift is
           ~1.5″/century versus ~22,000″/century for plain leapfrog.        */
        const double cbrt2 = 1.2599210498948732;
        const double w1 = 1.0/(2.0 - cbrt2), w0 = 1.0 - 2.0*w1;
        kdk_step(sim, w1*sim->dt);
        kdk_step(sim, w0*sim->dt);
        kdk_step(sim, w1*sim->dt);
    } else {
        kdk_step(sim, sim->dt);
    }
    sim->time_yr += sim->dt;
}

/* ── conserved quantities (Newtonian) ────────────────────────────────────── */
static void sim_energy(const Sim *sim, double *E, double *L){
    double K=0, U=0; Vec3 Lv={0,0,0};
    for(int i=0;i<N_BODIES;i++){
        K  += 0.5*sim->mass[i]*v3dot(sim->vel[i], sim->vel[i]);
        Lv  = v3add(Lv, v3scale(v3cross(sim->pos[i], sim->vel[i]), sim->mass[i]));
        for(int j=i+1;j<N_BODIES;j++)
            U -= G_AU*sim->mass[i]*sim->mass[j] / v3len(v3sub(sim->pos[j], sim->pos[i]));
    }
    *E = K+U; *L = v3len(Lv);
}

/* Osculating longitude of perihelion of body i (heliocentric), from the
   Laplace–Runge–Lenz / eccentricity vector  e = (v×h)/μ − r̂.  Returns the
   angle of e projected on the ecliptic, in radians.                       */
static double osc_varpi(const Sim *sim, int i, double *ecc_out){
    Vec3 r = v3sub(sim->pos[i], sim->pos[0]);
    Vec3 v = v3sub(sim->vel[i], sim->vel[0]);
    double mu = G_AU*(sim->mass[0]+sim->mass[i]);
    Vec3 h  = v3cross(r, v);
    Vec3 ev = v3sub(v3scale(v3cross(v, h), 1.0/mu), v3norm(r));
    if(ecc_out) *ecc_out = v3len(ev);
    return atan2(ev.y, ev.x);
}

/* ── simulation init ─────────────────────────────────────────────────────── */
static void sim_init(Sim *sim){
    memset(sim,0,sizeof(*sim));
    sim->dt         = DT_HOURS/(365.25*24.0);
    sim->integrator = INTEG_YOSHIDA4;
    sim->gr         = 0;

    for(int i=0;i<N_BODIES;i++){
        sim->mass[i] = BODIES[i].mass;
#ifdef HAVE_HORIZONS_J2000
        /* Osculating barycentric state at J2000.0; velocities AU/day → AU/yr */
        sim->pos[i] = (Vec3){ HORIZONS_J2000[i][0], HORIZONS_J2000[i][1], HORIZONS_J2000[i][2] };
        sim->vel[i] = (Vec3){ HORIZONS_J2000[i][3]*365.25, HORIZONS_J2000[i][4]*365.25,
                              HORIZONS_J2000[i][5]*365.25 };
#else
        /* The JPL table gives *mean* elements.  Used as osculating ones they
           reproduce Horizons to ~0.1 % for the inner planets but only to a
           few % for Saturn (Jupiter–Saturn great inequality).  Run
           tests/fetch_j2000_state.sh for exact initial conditions.       */
        elements_to_cartesian(&BODIES[i], &sim->pos[i], &sim->vel[i]);
#endif
    }
    /* Kick-drift-kick (velocity Verlet) keeps v synchronised with x, so the
       initial velocities are used as they are.  (An earlier version kicked
       them back by dt/2 — the initialisation for the *staggered* leapfrog —
       which biased every planet's eccentricity vector by ~a·dt/2.)        */
    compute_accel(sim->pos, sim->vel, sim->mass, sim->acc, sim->gr);
    sim_energy(sim, &sim->E0, &sim->L0);

    for(int i=0;i<N_BODIES;i++){
        if(BODIES[i].period_yr <= 0.0){
            sim->trail_cap[i] = 0;
        } else {
            double samples_per_yr = 365.25 / 3.0;      /* one every 3 days */
            int cap = (int)(BODIES[i].period_yr * samples_per_yr);
            if(cap < 50)        cap = 50;
            if(cap > TRAIL_LEN) cap = TRAIL_LEN;
            sim->trail_cap[i] = cap;
        }
    }
}

/* ── trail management ───────────────────────────────────────────────────── */
static void trail_push(Sim *sim, int i){
    int cap = sim->trail_cap[i];
    if(cap <= 0) return;
    int h = sim->trail_head[i];
    sim->trail[i][h] = sim->pos[i];
    sim->trail_head[i] = (h+1) % cap;
    if(sim->trail_count[i] < cap) sim->trail_count[i]++;
}

/* ── calendar date from simulation time (TT) ─────────────────────────────── */
static void sim_date(double time_yr, int *Y, int *Mo, int *D, int *h, int *mi){
    double jd = JD_J2000 + time_yr*365.25;
    double Z  = floor(jd + 0.5), F = jd + 0.5 - Z;
    double A  = Z;
    if(Z >= 2299161){ double al = floor((Z-1867216.25)/36524.25); A = Z+1+al-floor(al/4); }
    double B = A+1524, C = floor((B-122.1)/365.25), Dd = floor(365.25*C);
    double E = floor((B-Dd)/30.6001);
    *D  = (int)(B - Dd - floor(30.6001*E));
    *Mo = (int)(E < 14 ? E-1 : E-13);
    *Y  = (int)(*Mo > 2 ? C-4716 : C-4715);
    int secs = (int)floor(F*86400.0 + 0.5);
    *h = secs/3600; *mi = (secs%3600)/60;
}

/* ======================================================================
   Massless test particles (asteroid + Kuiper belts)
   Integrated with plain leapfrog at AST_SUB×dt under the massive bodies'
   current positions.  They feel the planets; the planets don't feel them.
   ====================================================================== */
typedef struct {
    Vec3 pos[N_AST], vel[N_AST], acc[N_AST];
    int  pending;             /* main steps since last asteroid step        */
    int  initialised;
} Belt;

static double frand(unsigned *s){           /* deterministic LCG, [0,1)   */
    *s = *s*1664525u + 1013904223u;
    return (*s >> 8) / 16777216.0;
}

static void belt_accel(const Sim *sim, Belt *b){
    for(int k=0;k<N_AST;k++){
        Vec3 a = {0,0,0};
        for(int j=0;j<N_BODIES;j++){
            Vec3 r = v3sub(sim->pos[j], b->pos[k]);
            double d2 = v3dot(r,r) + 1e-10, d3 = d2*sqrt(d2);
            a = v3add(a, v3scale(r, G_AU*sim->mass[j]/d3));
        }
        b->acc[k] = a;
    }
}

static void belt_init(const Sim *sim, Belt *b){
    unsigned seed = 20260130u;
    double mu = G_AU*sim->mass[0];
    for(int k=0;k<N_AST;k++){
        int kuiper = k >= N_AST_MAIN;
        double a   = kuiper ? 39.0 + 9.0*frand(&seed) : 2.1 + 1.2*frand(&seed);
        double e   = (kuiper ? 0.20 : 0.15) * sqrt(frand(&seed));       /* Rayleigh-ish */
        double inc = (kuiper ? 20.0 : 12.0) * frand(&seed) * M_PI/180.0;
        double Om  = 2*M_PI*frand(&seed), w = 2*M_PI*frand(&seed), M0 = 2*M_PI*frand(&seed);
        Vec3 p, v;
        kepler_to_cartesian(a, e, inc, Om, w, M0, mu, &p, &v);
        b->pos[k] = v3add(p, sim->pos[0]);                 /* heliocentric → barycentric */
        b->vel[k] = v3add(v, sim->vel[0]);
    }
    belt_accel(sim, b);
    b->pending = 0; b->initialised = 1;
}

/* Called once per main step; advances the belt every AST_SUB main steps. */
static void belt_tick(const Sim *sim, Belt *b){
    if(!b->initialised) return;
    if(++b->pending < AST_SUB) return;
    b->pending = 0;
    double h = AST_SUB*sim->dt, hh = 0.5*h;
    for(int k=0;k<N_AST;k++) b->vel[k] = v3add(b->vel[k], v3scale(b->acc[k], hh));
    for(int k=0;k<N_AST;k++) b->pos[k] = v3add(b->pos[k], v3scale(b->vel[k], h));
    belt_accel(sim, b);
    for(int k=0;k<N_AST;k++) b->vel[k] = v3add(b->vel[k], v3scale(b->acc[k], hh));
}

/* ======================================================================
   Camera / view
   ====================================================================== */
typedef struct {
    Mat3   rot;
    double scale;             /* logical px per AU (before perspective)     */
    double fov_dist;
    int    dragging;
    int    panning;
    double last_x, last_y;    /* last mouse position, logical px            */
    double pan_x, pan_y;      /* logical px                                 */
    int    selected;
} Camera;

static Camera cam_init(void){
    Camera c;
    c.rot      = mat3_identity();
    c.scale    = CAM_SCALE0;
    c.fov_dist = 8.0;
    c.dragging = 0;
    c.panning  = 0;
    c.last_x   = 0;
    c.last_y   = 0;
    c.pan_x    = 0.0;
    c.pan_y    = 0.0;
    c.selected = 3;
    Vec3 xax = {1,0,0};
    Mat3 tilt = mat3_rot_axis(xax, -0.4);
    c.rot = mat3_mul(&tilt, &c.rot);
    return c;
}

static void cam_set_view(Camera *c, int axis){
    Mat3 m = mat3_identity();
    switch(axis){
        case 0: m=mat3_rot_axis((Vec3){0,1,0}, -M_PI/2.0); break;
        case 1: m=mat3_rot_axis((Vec3){1,0,0},  M_PI/2.0); break;
        case 2: m=mat3_identity(); break;
    }
    c->rot = m;
}

/* Perspective weight for a view-space depth; project() applies it. */
static double cam_w(const Camera *cam, double vz){
    return cam->fov_dist / (cam->fov_dist + vz + 40.0);
}

/* World → logical screen coordinates. */
static void project(const Camera *cam, Vec3 world,
                    float *sx, float *sy, float *sz_out){
    Vec3 v = mat3_mul_vec(&cam->rot, world);
    double w = cam_w(cam, v.z);
    *sx = (float)(SIM_W/2.0 + cam->pan_x + v.x * cam->scale * w);
    *sy = (float)(WIN_H/2.0 + cam->pan_y - v.y * cam->scale * w);
    if(sz_out) *sz_out = (float)v.z;
}

/* ======================================================================
   Device transform
   The scene is laid out in logical WIN_W×WIN_H; the back-buffer is the
   display's physical pixel size.  g_S scales, g_ox/g_oy letterbox when the
   aspect differs (native window), so every drawn pixel is a device pixel:
   no resampling in either direction, hence crisp at any size and DPR.
   ====================================================================== */
static double g_S = 1.0, g_ox = 0.0, g_oy = 0.0;
static int    g_W = WIN_W, g_H = WIN_H;      /* back-buffer, device px      */

static inline float DX(double x){ return (float)(g_ox + x*g_S); }
static inline float DY(double y){ return (float)(g_oy + y*g_S); }
static inline int   DL(double v){ int t=(int)lround(v*g_S); return t<1?1:t; } /* a length */
static inline int   stroke_px(void){ return DL(1.0); }

/* ── primitives (logical in, device out) ─────────────────────────────── */
static void set_color(SDL_Renderer *rend, Uint8 r, Uint8 g, Uint8 b, Uint8 a){
    SDL_SetRenderDrawColor(rend, r, g, b, a);
}

/* 1-logical-px dot → stroke_px square */
static void dot(SDL_Renderer *rend, double x, double y){
    int t = stroke_px();
    SDL_Rect rc = { (int)lround(DX(x)) - t/2, (int)lround(DY(y)) - t/2, t, t };
    SDL_RenderFillRect(rend, &rc);
}

/* Line with stroke weight round(g_S): thickness copies offset across the
   minor axis.  Drawn in device space so it is continuous at any scale.   */
static void line(SDL_Renderer *rend, double x0, double y0, double x1, double y1){
    float X0=DX(x0), Y0=DY(y0), X1=DX(x1), Y1=DY(y1);
    int t = stroke_px();
    int steep = fabsf(Y1-Y0) > fabsf(X1-X0);
    for(int k=0;k<t;k++){
        float o = (float)(k - t/2);
        if(steep) SDL_RenderDrawLineF(rend, X0+o, Y0, X1+o, Y1);
        else      SDL_RenderDrawLineF(rend, X0, Y0+o, X1, Y1+o);
    }
}

/* Device-space polyline, n points already converted with DX/DY.
   One SDL_RenderDrawLinesF call per thickness copy: SDL paints interior
   joints once, so alpha-blended curves don't get a bright dot per vertex
   the way per-segment calls do.                                           */
static void polyline_dev(SDL_Renderer *rend, const SDL_FPoint *p, int n){
    static SDL_FPoint tmp[TRAIL_LEN + 2];
    if(n < 2) return;
    int t = stroke_px();
    for(int k=0;k<t;k++){
        float o = (float)(k - t/2);
        for(int i=0;i<n;i++){
            int a = i-8 < 0 ? 0 : i-8, b = i+8 > n-1 ? n-1 : i+8;
            int steep = fabsf(p[b].y-p[a].y) > fabsf(p[b].x-p[a].x);
            tmp[i].x = p[i].x + (steep ? o : 0);
            tmp[i].y = p[i].y + (steep ? 0 : o);
        }
        SDL_RenderDrawLinesF(rend, tmp, n);
    }
}

static void fill_circle(SDL_Renderer *rend, double cx, double cy, double rad){
    int R = DL(rad); int CX=(int)lround(DX(cx)), CY=(int)lround(DY(cy));
    for(int dy=-R; dy<=R; dy++){
        int dx = (int)sqrt((double)(R*R - dy*dy));
        SDL_RenderDrawLine(rend, CX-dx, CY+dy, CX+dx, CY+dy);
    }
}

static void ring(SDL_Renderer *rend, double cx, double cy, double rad){
    double px=cx+rad, py=cy;
    for(int ang=6; ang<=360; ang+=6){
        double a=ang*M_PI/180.0, x=cx+rad*cos(a), y=cy+rad*sin(a);
        line(rend, px, py, x, y); px=x; py=y;
    }
}

/* ── text: fonts are opened at size×g_S, so glyphs are rendered by FreeType
   at device resolution and blitted 1:1.  Positions are logical.          */
static int text_width(TTF_Font *font, const char *text){   /* logical px */
    int w=0, h=0;
    TTF_SizeUTF8(font, text, &w, &h);
    return (int)ceil(w / g_S);
}

static void draw_text(SDL_Renderer *rend, TTF_Font *font,
                      const char *text, double x, double y,
                      Uint8 r, Uint8 g, Uint8 b){
    SDL_Color col = {r,g,b,255};
    SDL_Surface *surf = TTF_RenderUTF8_Blended(font, text, col);
    if(!surf) return;
    SDL_Texture *tex = SDL_CreateTextureFromSurface(rend, surf);
    if(tex){
        SDL_Rect dst = {(int)lround(DX(x)), (int)lround(DY(y)), surf->w, surf->h};
        SDL_RenderCopy(rend, tex, NULL, &dst);
        SDL_DestroyTexture(tex);
    }
    SDL_FreeSurface(surf);
}

/* ======================================================================
   Application state — file-scope so loop_iter() works from both the native
   while-loop and Emscripten's callback-driven main loop.
   ====================================================================== */
typedef struct {            /* linear fit of unwrapped ϖ(t) → secular rate */
    int    body;
    double t0, w_prev, w_unw, w0;
    double n, sx, sy, sxx, sxy, syy;
} Precession;

static struct {
    SDL_Window   *win;
    SDL_Renderer *rend;
    TTF_Font     *font_lg, *font_sm;
    char          font_path[512];
    int           font_px_lg, font_px_sm;
    Camera        cam;
    int           trail_tick;
    int           steps_per_frame;
    int           running;
    Uint32        last_frame;
    double        fps;
    double        mouse_k;      /* device px per window unit (HiDPI)       */

    int paused, follow, show_help, show_belt;
    double E, L;                /* current energy / |L| for the readout     */
    Precession prec;
#ifndef __EMSCRIPTEN__
    FILE  *ffmpeg_pipe;
    Uint8 *pixel_buf;
    int    recording, rec_frames, rec_w, rec_h;
#endif
} app;

static Sim  g_sim;          /* ~6 MB of trail buffers — keep it static     */
static Belt g_belt;

/* ── perihelion-drift tracker ────────────────────────────────────────── */
static void prec_reset(Precession *p, const Sim *sim, int body){
    memset(p, 0, sizeof *p);
    p->body = body; p->t0 = sim->time_yr;
    p->w0 = p->w_prev = osc_varpi(sim, body, NULL);
}
static void prec_sample(Precession *p, const Sim *sim){
    double w = osc_varpi(sim, p->body, NULL);
    double d = w - p->w_prev;
    if(d >  M_PI) d -= 2*M_PI;
    if(d < -M_PI) d += 2*M_PI;
    p->w_unw += d; p->w_prev = w;
    double x = sim->time_yr - p->t0, y = p->w_unw;
    p->n++; p->sx+=x; p->sy+=y; p->sxx+=x*x; p->sxy+=x*y; p->syy+=y*y;
}
/* Returns 1 and fills rate ± se in arcsec/century once enough data exists.
   Short-period perturbations (Venus, Jupiter) dominate the slope for the
   first few years, so nothing is reported before 5 yr of simulated time. */
static int prec_rate(const Precession *p, double *rate, double *se, double *span){
    double n=p->n; if(n < 100) return 0;
    if(fabs(p->sx/n)*2.0 < 5.0) return 0;
    double det = n*p->sxx - p->sx*p->sx;
    if(fabs(det) < 1e-12) return 0;
    double b  = (n*p->sxy - p->sx*p->sy)/det;
    double ss = p->syy - p->sy*p->sy/n - b*b*(p->sxx - p->sx*p->sx/n);
    if(ss < 0) ss = 0;
    double sxx_c = p->sxx - p->sx*p->sx/n;
    *rate = b * 206264.806 * 100.0;
    *se   = sqrt(ss/(n-2)/sxx_c) * 206264.806 * 100.0;
    *span = fabs(p->sx/n)*2.0;   /* ≈ time span covered */
    return 1;
}

/* ======================================================================
   Scene rendering
   ====================================================================== */
static void draw_orbit(SDL_Renderer *rend, const Camera *cam, const BodyDef *bd){
    static SDL_FPoint pts[361];
    if(bd->a == 0.0) return;
    double d = M_PI/180.0;
    double inc=bd->inc*d, Om=bd->omega*d, w=bd->w*d;
    double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);
    double Qxx = cO*cw - sO*sw*ci,  Qxy = -cO*sw - sO*cw*ci;
    double Qyx = sO*cw + cO*sw*ci,  Qyy = -sO*sw + cO*cw*ci;
    double Qzx = sw*si,              Qzy =  cw*si;
    int n=0;
    for(int step=0; step<=360; step++){
        double nu   = 2.0*M_PI*step/360.0;
        double dist = bd->a*(1.0 - bd->ecc*bd->ecc)/(1.0 + bd->ecc*cos(nu));
        double ox = dist*cos(nu), oy = dist*sin(nu);
        Vec3 p = {Qxx*ox+Qxy*oy, Qyx*ox+Qyy*oy, Qzx*ox+Qzy*oy};
        float sx,sy; project(cam, v3add(p, g_sim.pos[0]), &sx, &sy, NULL);
        pts[n].x = DX(sx); pts[n].y = DY(sy); n++;
    }
    polyline_dev(rend, pts, n);
}

static void draw_trail(SDL_Renderer *rend, const Camera *cam, const Sim *sim, int i){
    static SDL_FPoint pts[TRAIL_LEN + 1];
    int cnt = sim->trail_count[i], head = sim->trail_head[i], cap = sim->trail_cap[i];
    if(cap <= 0 || cnt < 2) return;
    /* Alpha fades along the trail.  Draw it as ~16 constant-alpha runs that
       share end points, so joints are painted once (see polyline_dev).  */
    int runs = cnt < 16 ? 1 : 16, per = (cnt + runs - 1)/runs;
    for(int r=0; r<runs; r++){
        int k0 = r*per, k1 = (r+1)*per; if(k1 > cnt-1) k1 = cnt-1;
        if(k1 <= k0) break;
        int n=0;
        for(int k=k0; k<=k1; k++){
            int idx = ((head - cnt + k) + cap) % cap;
            float sx,sy; project(cam, sim->trail[i][idx], &sx, &sy, NULL);
            pts[n].x = DX(sx); pts[n].y = DY(sy); n++;
        }
        Uint8 alpha = (Uint8)(180.0 * (k0 + k1) / (2.0*cnt));
        SDL_SetRenderDrawColor(rend, BODIES[i].r, BODIES[i].g, BODIES[i].b, alpha);
        polyline_dev(rend, pts, n);
    }
}

static void draw_belt(SDL_Renderer *rend, const Camera *cam){
    if(!app.show_belt || !g_belt.initialised) return;
    int t = stroke_px();
    for(int k=0;k<N_AST;k++){
        float sx,sy; project(cam, g_belt.pos[k], &sx, &sy, NULL);
        if(sx < 0 || sx >= SIM_W || sy < 0 || sy >= WIN_H) continue;
        SDL_Rect rc = { (int)lround(DX(sx)) - t/2, (int)lround(DY(sy)) - t/2, t, t };
        if(k < N_AST_MAIN) SDL_SetRenderDrawColor(rend, 170,150,120, 150);
        else               SDL_SetRenderDrawColor(rend, 120,150,190, 130);
        SDL_RenderFillRect(rend, &rc);
    }
}

/* Labels: greedy placement right / below / above of the body; a label that
   would overlap one already placed is dropped — except for the Sun and
   the selected body, which are always shown.                             */
typedef struct { double x, y, w, h; } Box;
static int box_hit(Box a, Box b){
    return !(a.x+a.w < b.x || b.x+b.w < a.x || a.y+a.h < b.y || b.y+b.h < a.y);
}

static void draw_bodies(SDL_Renderer *rend, const Camera *cam, const Sim *sim){
    Box placed[N_BODIES]; int nplaced = 0;
    /* Sun first so inner planets yield to it; then far-to-near by depth so
       the nearest bodies win label slots. */
    int order[N_BODIES]; float depth[N_BODIES];
    for(int i=0;i<N_BODIES;i++){ float sx,sy; project(cam, sim->pos[i], &sx,&sy,&depth[i]); order[i]=i; }
    for(int a=1;a<N_BODIES;a++) for(int b=a+1;b<N_BODIES;b++)
        if(depth[order[b]] > depth[order[a]]){ int t=order[a]; order[a]=order[b]; order[b]=t; }

    for(int oi=0; oi<N_BODIES; oi++){
        int i = order[oi];
        float sx,sy; project(cam, sim->pos[i], &sx, &sy, NULL);
        double rad = BODIES[i].draw_r;

        if(i==0){   /* corona */
            for(int glow=(int)rad+8; glow>(int)rad; glow--){
                SDL_SetRenderDrawColor(rend,255,200,50,(Uint8)(30*(glow-rad)));
                fill_circle(rend, sx, sy, glow);
            }
        }
        set_color(rend, BODIES[i].r, BODIES[i].g, BODIES[i].b, 255);
        fill_circle(rend, sx, sy, rad);

        if(i == cam->selected){
            SDL_SetRenderDrawColor(rend,255,255,255,160);
            ring(rend, sx, sy, rad+4);
        }

        /* label */
        int tw = text_width(app.font_sm, BODIES[i].name), th = 14;
        Box cands[3] = {
            { sx+rad+4,     sy-7,          tw, th },
            { sx-tw/2.0,    sy+rad+4,      tw, th },
            { sx-tw/2.0,    sy-rad-4-th,   tw, th },
        };
        int chosen = -1;
        for(int c=0;c<3 && chosen<0;c++){
            int ok = 1;
            for(int p=0;p<nplaced;p++) if(box_hit(cands[c], placed[p])){ ok=0; break; }
            if(ok) chosen = c;
        }
        if(chosen < 0 && (i==0 || i==cam->selected)) chosen = 0;
        if(chosen >= 0){
            placed[nplaced++] = cands[chosen];
            draw_text(rend, app.font_sm, BODIES[i].name, cands[chosen].x, cands[chosen].y,
                      BODIES[i].r, BODIES[i].g, BODIES[i].b);
        }
    }
}

static void draw_axes(SDL_Renderer *rend, const Camera *cam){
    double ax0=50, ay0=WIN_H-50, len=30.0;
    Vec3 ex={1,0,0},ey={0,1,0},ez={0,0,1};
    Vec3 rx=mat3_mul_vec(&cam->rot,ex), ry=mat3_mul_vec(&cam->rot,ey), rz=mat3_mul_vec(&cam->rot,ez);
    SDL_SetRenderDrawColor(rend,200,50,50,200);
    line(rend,ax0,ay0, ax0+len*rx.x, ay0-len*rx.y);
    draw_text(rend,app.font_sm,"X",ax0+len*rx.x+2,ay0-len*rx.y,200,50,50);
    SDL_SetRenderDrawColor(rend,50,200,50,200);
    line(rend,ax0,ay0, ax0+len*ry.x, ay0-len*ry.y);
    draw_text(rend,app.font_sm,"Y",ax0+len*ry.x+2,ay0-len*ry.y,50,200,50);
    SDL_SetRenderDrawColor(rend,50,100,220,200);
    line(rend,ax0,ay0, ax0+len*rz.x, ay0-len*rz.y);
    draw_text(rend,app.font_sm,"Z",ax0+len*rz.x+2,ay0-len*rz.y,50,100,220);
}

/* ── info panel ───────────────────────────────────────────────────────── */
static void draw_panel(SDL_Renderer *rend, const Sim *sim, const Camera *cam){
    TTF_Font *font_lg = app.font_lg, *font_sm = app.font_sm;
    int sel = cam->selected;
    double px0 = SIM_W + 8;
    char buf[160];

    set_color(rend,15,15,25,255);
    SDL_Rect bg = {(int)lround(DX(SIM_W)), (int)lround(DY(0)),
                   (int)lround(PANEL_W*g_S)+1, (int)lround(WIN_H*g_S)+1};
    SDL_RenderFillRect(rend, &bg);
    set_color(rend,60,60,80,255);
    line(rend, SIM_W, 0, SIM_W, WIN_H);

    draw_text(rend, font_lg, "Solar System", px0, 8, 255,220,50);

    /* ── header: date, speed, conservation, integrator, precession ─── */
    int Y,Mo,D,h,mi; sim_date(sim->time_yr, &Y,&Mo,&D,&h,&mi);
    snprintf(buf, sizeof buf, "%04d-%02d-%02d %02d:%02d TT   (T = %+.2f yr)", Y,Mo,D,h,mi, sim->time_yr);
    draw_text(rend, font_sm, buf, px0, 36, 180,180,200);

    double days_per_s = app.steps_per_frame * fabs(sim->dt) * 365.25 * app.fps;
    if(app.paused)
        snprintf(buf, sizeof buf, "PAUSED   (%d steps/frame)", app.steps_per_frame);
    else if(days_per_s >= 365.25)
        snprintf(buf, sizeof buf, "Speed: %s%.1f yr/s  (%d steps/frame)",
                 sim->dt<0?"\xe2\x97\x80 ":"", days_per_s/365.25, app.steps_per_frame);
    else
        snprintf(buf, sizeof buf, "Speed: %s%.1f d/s  (%d steps/frame)",
                 sim->dt<0?"\xe2\x97\x80 ":"", days_per_s, app.steps_per_frame);
    draw_text(rend, font_sm, buf, px0, 52, 180,180,200);

    snprintf(buf, sizeof buf, "\xce\x94""E/E\xe2\x82\x80 = %+.1e    \xce\x94""L/L\xe2\x82\x80 = %+.1e",
             (app.E - sim->E0)/fabs(sim->E0), (app.L - sim->L0)/sim->L0);
    draw_text(rend, font_sm, buf, px0, 68, sim->gr ? 140:160, 160, sim->gr ? 140:120);

    snprintf(buf, sizeof buf, "%s, dt = %g h   GR (1PN): %s%s",
             sim->integrator==INTEG_YOSHIDA4 ? "Yoshida-4" : "Leapfrog-2",
             DT_HOURS, sim->gr ? "on" : "off",
             sim->gr ? "  (E is Newtonian)" : "");
    draw_text(rend, font_sm, buf, px0, 84, 140,140,160);

    {
        double rate, se, span;
        if(prec_rate(&app.prec, &rate, &se, &span))
            snprintf(buf, sizeof buf, "Mercury \xcf\x96\xcc\x87 = %+.1f \xc2\xb1 %.1f \xe2\x80\xb3/cy  (%.1f yr)",
                     rate, se, span);
        else
            snprintf(buf, sizeof buf, "Mercury \xcf\x96\xcc\x87: collecting\xe2\x80\xa6 (%.1f of 5 yr)",
                     fabs(sim->time_yr - app.prec.t0));
        draw_text(rend, font_sm, buf, px0, 100, 180,160,100);
    }

    set_color(rend,50,50,70,255);
    line(rend, SIM_W+4, 120, WIN_W-4, 120);

    /* ── body list ─────────────────────────────────────────────────── */
    double y0 = 128;
    for(int i=0;i<N_BODIES;i++){
        const BodyDef *bd = &BODIES[i];
        if(i == sel){
            SDL_Rect hl = {(int)lround(DX(SIM_W)), (int)lround(DY(y0-2)),
                           (int)lround(PANEL_W*g_S), (int)lround(22*g_S)};
            set_color(rend,30,30,55,255);
            SDL_RenderFillRect(rend, &hl);
        }
        set_color(rend, bd->r, bd->g, bd->b, 255);
        fill_circle(rend, px0+6, y0+8, 5);
        draw_text(rend, font_sm, bd->name, px0+16, y0,
                  (i==sel)?255:200, (i==sel)?255:200, (i==sel)?255:200);
        y0 += 22;
    }

    set_color(rend,50,50,70,255);
    line(rend, SIM_W+4, y0+4, WIN_W-4, y0+4);

    /* ── selected body ─────────────────────────────────────────────── */
    if(sel >= 0 && sel < N_BODIES){
        const BodyDef *bd = &BODIES[sel];
        double yd = y0 + 14;
        snprintf(buf, sizeof buf, "%s%s", bd->name, app.follow ? "   (following)" : "");
        draw_text(rend, font_lg, buf, px0, yd, bd->r, bd->g, bd->b);
        yd += 26;

        Vec3 r = v3sub(sim->pos[sel], sim->pos[0]);
        double dist_AU = v3len(r);
        Vec3 vrel  = v3sub(sim->vel[sel], sim->vel[0]);
        double v_km_s = v3len(vrel) * 4.740470;          /* AU/yr → km/s   */
        double e_osc = 0.0; if(sel>0) osc_varpi(sim, sel, &e_osc);

        char v_sma[32], v_ecc[32], v_inc[32], v_per[32];
        char v_dist[32], v_spd[32], v_mspd[32], v_rad[32], v_mass[32], v_eosc[32];
        snprintf(v_sma,  sizeof v_sma,  "%.4f AU",   bd->a);
        snprintf(v_ecc,  sizeof v_ecc,  "%.5f",      bd->ecc);
        snprintf(v_inc,  sizeof v_inc,  "%.3f\xc2\xb0", bd->inc);
        if(bd->period_yr > 0) snprintf(v_per, sizeof v_per, "%.3f yr", bd->period_yr);
        else                  snprintf(v_per, sizeof v_per, "\xe2\x80\x93");
        snprintf(v_dist, sizeof v_dist, "%.4f AU",   dist_AU);
        snprintf(v_spd,  sizeof v_spd,  "%.2f km/s", v_km_s);
        snprintf(v_mspd, sizeof v_mspd, "%.2f km/s", bd->v_orb_km_s);
        if(sel>0) snprintf(v_eosc, sizeof v_eosc, "%.5f", e_osc); else v_eosc[0]=0;
        snprintf(v_rad,  sizeof v_rad,  "%.0f km",   bd->radius_km);
        snprintf(v_mass, sizeof v_mass, "%.3e M\xe2\x8a\x99", bd->mass);

        struct { const char *lbl; const char *val; int section; } rows[] = {
            {"Orbital Elements (J2000)", NULL,   1},
            {"  Semi-major axis",  v_sma,   0},
            {"  Eccentricity",     v_ecc,   0},
            {"  Inclination",      v_inc,   0},
            {"  Period",           v_per,   0},
            {"Current State",      NULL,    1},
            {"  Dist. from Sun",   v_dist,  0},
            {"  Speed (now)",      v_spd,   0},
            {"  Speed (mean)",     v_mspd,  0},
            {"  Osc. eccentricity",v_eosc,  0},
            {"Physical",           NULL,    1},
            {"  Radius",           v_rad,   0},
            {"  Mass",             v_mass,  0},
        };
        int val_x = 128;
        for(int k=0; k<(int)(sizeof(rows)/sizeof(rows[0])); k++){
            if(rows[k].section){
                draw_text(rend, font_sm, rows[k].lbl, px0, yd, 180,180,100);
                yd += 18;
            } else {
                if(!rows[k].val || !rows[k].val[0]) continue;
                draw_text(rend, font_sm, rows[k].lbl, px0, yd, 140,140,160);
                draw_text(rend, font_sm, rows[k].val, px0+val_x, yd, 220,220,240);
                yd += 16;
            }
        }
    }

    /* ── controls hint ─────────────────────────────────────────────── */
    double yh = WIN_H - 58;
    set_color(rend,50,50,70,255);
    line(rend, SIM_W+4, yh-4, WIN_W-4, yh-4);
    draw_text(rend, font_sm, "H  help \xe2\x80\x93 all keys",                   px0, yh,    160,160,100);
    draw_text(rend, font_sm, "Drag rotate \xc2\xb7 Shift+drag pan \xc2\xb7 Scroll zoom \xc2\xb7 Click select", px0, yh+14, 120,120,140);
    draw_text(rend, font_sm, "Space pause \xc2\xb7 T reverse \xc2\xb7 +/- speed \xc2\xb7 F follow",    px0, yh+28, 120,120,140);
    draw_text(rend, font_sm, "G relativity \xc2\xb7 I integrator \xc2\xb7 A asteroids \xc2\xb7 1-9 body", px0, yh+42, 120,120,140);
}

/* ── help overlay ─────────────────────────────────────────────────────── */
static void draw_help(SDL_Renderer *rend){
    static const char *rows[][2] = {
        {"Drag",        "rotate view (arcball)"},
        {"Shift+drag",  "pan"},
        {"Scroll",      "zoom in / out"},
        {"Click body",  "select it"},
        {"1 \xe2\x80\x93 9",   "select Sun \xe2\x80\xa6 Neptune"},
        {"F",           "follow selected body (centre + zoom to its orbit)"},
        {"Z / X / Y",   "top-down (XY) / edge-on (YZ) / edge-on (XZ) view"},
        {"R",           "reset view, zoom and pan"},
        {"Space",       "pause / resume"},
        {"T",           "reverse time (the integrator is time-symmetric)"},
        {"+ / \xe2\x88\x92",   "faster / slower (integration steps per frame)"},
        {"G",           "general relativity: 1PN Sun\xe2\x80\x93planet term on/off"},
        {"I",           "integrator: Yoshida-4 \xe2\x87\x84 Leapfrog-2"},
        {"A",           "asteroid + Kuiper belt test particles on/off"},
#ifndef __EMSCRIPTEN__
        {"V",           "start / stop MP4 capture (ffmpeg)"},
#endif
        {"H",           "this help"},
        {"Q / Esc",     "quit"},
    };
    int n = (int)(sizeof rows / sizeof rows[0]);
    double pw = 560, ph = 70 + n*20 + 56;
    double px = (SIM_W - pw)/2, py = (WIN_H - ph)/2;

    SDL_SetRenderDrawColor(rend, 0,0,0, 170);
    SDL_Rect dim = {(int)lround(DX(0)), (int)lround(DY(0)), (int)lround(SIM_W*g_S), (int)lround(WIN_H*g_S)};
    SDL_RenderFillRect(rend, &dim);
    SDL_SetRenderDrawColor(rend, 12,14,30, 235);
    SDL_Rect panel = {(int)lround(DX(px)), (int)lround(DY(py)), (int)lround(pw*g_S), (int)lround(ph*g_S)};
    SDL_RenderFillRect(rend, &panel);
    SDL_SetRenderDrawColor(rend, 80,110,200, 255);
    line(rend, px,py, px+pw,py); line(rend, px+pw,py, px+pw,py+ph);
    line(rend, px+pw,py+ph, px,py+ph); line(rend, px,py+ph, px,py);

    draw_text(rend, app.font_lg, "Key bindings", px+18, py+14, 100,200,255);
    double y = py+50;
    for(int i=0;i<n;i++, y+=20){
        draw_text(rend, app.font_sm, rows[i][0], px+18,  y, 120,200,255);
        draw_text(rend, app.font_sm, rows[i][1], px+130, y, 220,220,240);
    }
    draw_text(rend, app.font_sm,
        "Readouts: \xce\x94""E, \xce\x94L are drifts of the conserved quantities since t=0;",
        px+18, y+10, 150,150,170);
    draw_text(rend, app.font_sm,
        "\xcf\x96\xcc\x87 is Mercury's secular perihelion rate (GR adds 43\xe2\x80\xb3/century).",
        px+18, y+26, 150,150,170);
}

/* ======================================================================
   Resolution management
   ====================================================================== */
static int open_fonts(double S){
    int lg = (int)lround(15*S), sm = (int)lround(12*S);
    if(lg < 6) lg = 6;
    if(sm < 5) sm = 5;
    if(app.font_lg && lg == app.font_px_lg && sm == app.font_px_sm) return 1;
    TTF_Font *flg = TTF_OpenFont(app.font_path, lg);
    TTF_Font *fsm = TTF_OpenFont(app.font_path, sm);
    if(!flg || !fsm){ if(flg) TTF_CloseFont(flg); if(fsm) TTF_CloseFont(fsm); return 0; }
    if(app.font_lg) TTF_CloseFont(app.font_lg);
    if(app.font_sm) TTF_CloseFont(app.font_sm);
    app.font_lg = flg; app.font_sm = fsm;
    app.font_px_lg = lg; app.font_px_sm = sm;
    return 1;
}

/* Recompute the logical→device mapping from the actual drawable size. */
static void apply_resolution(void){
    int W, H, ww, wh;
    SDL_GetRendererOutputSize(app.rend, &W, &H);
    SDL_GetWindowSize(app.win, &ww, &wh);
    g_W = W; g_H = H;
    app.mouse_k = ww > 0 ? (double)W/ww : 1.0;
    double sx = (double)W/WIN_W, sy = (double)H/WIN_H;
    g_S  = sx < sy ? sx : sy;
    g_ox = (W - WIN_W*g_S)/2.0;
    g_oy = (H - WIN_H*g_S)/2.0;
    open_fonts(g_S);
#ifdef __EMSCRIPTEN__
    EM_ASM({ if (window._fitCss) window._fitCss($0, $1); }, W, H);
#endif
}

#ifdef __EMSCRIPTEN__
static int g_last_devw = 0;
/* The page measures its container in device pixels (debounced); when it
   changes, resize the canvas back-buffer to match at the 16:9 layout aspect. */
static void poll_resize(void){
    int dw = EM_ASM_INT({ return (window._devW | 0); });
    if(dw <= 0 || dw == g_last_devw) return;
    g_last_devw = dw;
    if(dw < 640) dw = 640;
    if(dw > 3840) dw = 3840;
    int W = dw & ~1, H = (int)lround((double)W*WIN_H/WIN_W) & ~1;
    SDL_SetWindowSize(app.win, W, H);
    apply_resolution();
}
#endif

/* ── window units → logical coordinates ───────────────────────────────── */
static void to_logical(int wx, int wy, double *lx, double *ly){
    *lx = (wx*app.mouse_k - g_ox)/g_S;
    *ly = (wy*app.mouse_k - g_oy)/g_S;
}

/* ── follow mode: re-fit zoom to the selected body's orbit ─────────────── */
static void follow_fit(void){
    int s = app.cam.selected;
    if(s <= 0){ app.cam.scale = CAM_SCALE0; return; }
    Vec3 v = mat3_mul_vec(&app.cam.rot, g_sim.pos[s]);
    double w = cam_w(&app.cam, v.z);
    double r_max = BODIES[s].a*(1.0+BODIES[s].ecc);
    double half  = 0.42 * (SIM_W < WIN_H ? SIM_W : WIN_H);
    app.cam.scale = half / (r_max * w);
    if(app.cam.scale < 3.0)    app.cam.scale = 3.0;
    if(app.cam.scale > 8000.0) app.cam.scale = 8000.0;
}
static void follow_center(void){
    int s = app.cam.selected;
    Vec3 v = mat3_mul_vec(&app.cam.rot, g_sim.pos[s]);
    double w = cam_w(&app.cam, v.z);
    app.cam.pan_x = -v.x*app.cam.scale*w;
    app.cam.pan_y =  v.y*app.cam.scale*w;
}
static void select_body(int i){
    if(i < 0 || i >= N_BODIES) return;
    app.cam.selected = i;
    if(app.follow) follow_fit();
}

/* ── video capture (native) ──────────────────────────────────────────── */
#ifndef __EMSCRIPTEN__
static void video_toggle(void){
    if(!app.recording){
        time_t t = time(NULL);
        struct tm *tm_ = localtime(&t);
        char fname[64], cmd[256];
        strftime(fname, sizeof fname, "solar_%Y%m%d_%H%M%S.mp4", tm_);
        app.rec_w = g_W & ~1; app.rec_h = g_H & ~1;
        snprintf(cmd, sizeof cmd,
            "ffmpeg -y -f rawvideo -pixel_format bgra -video_size %dx%d -framerate 60 "
            "-i pipe:0 -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p %s 2>/dev/null",
            app.rec_w, app.rec_h, fname);
        app.ffmpeg_pipe = popen(cmd, "w");
        if(app.ffmpeg_pipe){
            app.pixel_buf = malloc((size_t)app.rec_w*app.rec_h*4);
            app.recording = 1; app.rec_frames = 0;
            SDL_Log("Recording started: %s (%dx%d)", fname, app.rec_w, app.rec_h);
        }
    } else {
        app.recording = 0;
        free(app.pixel_buf); app.pixel_buf = NULL;
        pclose(app.ffmpeg_pipe); app.ffmpeg_pipe = NULL;
        SDL_Log("Recording stopped (%d frames)", app.rec_frames);
    }
}
static void video_frame(void){
    if(!(app.recording && app.ffmpeg_pipe && app.pixel_buf)) return;
    if((g_W & ~1) != app.rec_w || (g_H & ~1) != app.rec_h) return;  /* resized */
    SDL_Rect rc = {0,0,app.rec_w,app.rec_h};
    SDL_RenderReadPixels(app.rend, &rc, SDL_PIXELFORMAT_BGRA32, app.pixel_buf, app.rec_w*4);
    fwrite(app.pixel_buf, 4, (size_t)app.rec_w*app.rec_h, app.ffmpeg_pipe);
    app.rec_frames++;
}

/* Headless debugging aid:
     SOLAR_SNAPSHOT=frame.bmp   render N frames (SOLAR_SNAPSHOT_FRAMES, default
                                30), save the frame, exit
     SOLAR_SNAPSHOT_KEYS=z5     keys to feed first (e.g. z top view, 5 Mars)
   Combine with SDL_VIDEODRIVER=dummy for a fully headless render.        */
static void maybe_snapshot(void){
    static int frames = 0, checked = 0, want = 30;
    static const char *path = NULL;
    if(!checked){
        checked = 1;
        path = getenv("SOLAR_SNAPSHOT");
        if(path){
            const char *nf = getenv("SOLAR_SNAPSHOT_FRAMES");
            if(nf && atoi(nf) > 0) want = atoi(nf);
            const char *keys = getenv("SOLAR_SNAPSHOT_KEYS");
            for(; keys && *keys; keys++){
                SDL_Event e; memset(&e, 0, sizeof e);
                e.type = SDL_KEYDOWN; e.key.keysym.sym = (SDL_Keycode)*keys;
                SDL_PushEvent(&e);
            }
        }
    }
    if(!path || ++frames < want) return;
    SDL_Surface *sf = SDL_CreateRGBSurfaceWithFormat(0, g_W, g_H, 32, SDL_PIXELFORMAT_ARGB8888);
    if(sf && SDL_RenderReadPixels(app.rend, NULL, SDL_PIXELFORMAT_ARGB8888, sf->pixels, sf->pitch) == 0)
        SDL_SaveBMP(sf, path);
    if(sf) SDL_FreeSurface(sf);
    app.running = 0;
}
#endif

static void request_quit(void){
#ifdef __EMSCRIPTEN__
    emscripten_cancel_main_loop();
#else
    app.running = 0;
#endif
}

/* ── keyboard ─────────────────────────────────────────────────────────── */
static void handle_key(SDL_Keycode k){
    Camera *cam = &app.cam;
    switch(k){
    case SDLK_ESCAPE: case SDLK_q: request_quit(); break;
    case SDLK_EQUALS: case SDLK_PLUS: case SDLK_KP_PLUS:
        app.steps_per_frame = app.steps_per_frame < 10 ? app.steps_per_frame+1
                            : (int)(app.steps_per_frame*1.25);
        if(app.steps_per_frame > 6000) app.steps_per_frame = 6000;
        break;
    case SDLK_MINUS: case SDLK_KP_MINUS:
        app.steps_per_frame = app.steps_per_frame <= 10 ? app.steps_per_frame-1
                            : (int)(app.steps_per_frame/1.25);
        if(app.steps_per_frame < 1) app.steps_per_frame = 1;
        break;
    case SDLK_SPACE: app.paused = !app.paused; break;
    case SDLK_t:     g_sim.dt = -g_sim.dt; break;
    case SDLK_r:
        app.follow = 0; app.cam = cam_init();
        break;
    case SDLK_x: cam_set_view(cam, 0); break;
    case SDLK_y: cam_set_view(cam, 1); break;
    case SDLK_z: cam_set_view(cam, 2); break;
    case SDLK_f:
        app.follow = !app.follow;
        if(app.follow) follow_fit();
        break;
    case SDLK_g:
        g_sim.gr = !g_sim.gr;
        compute_accel(g_sim.pos, g_sim.vel, g_sim.mass, g_sim.acc, g_sim.gr);
        prec_reset(&app.prec, &g_sim, 1);
        break;
    case SDLK_i:
        g_sim.integrator = g_sim.integrator==INTEG_YOSHIDA4 ? INTEG_LEAPFROG : INTEG_YOSHIDA4;
        prec_reset(&app.prec, &g_sim, 1);
        break;
    case SDLK_a:
        app.show_belt = !app.show_belt;
        if(app.show_belt && !g_belt.initialised) belt_init(&g_sim, &g_belt);
        break;
    case SDLK_h: app.show_help = !app.show_help; break;
    case SDLK_1: case SDLK_2: case SDLK_3: case SDLK_4: case SDLK_5:
    case SDLK_6: case SDLK_7: case SDLK_8: case SDLK_9:
        select_body((int)(k - SDLK_1)); break;
#ifndef __EMSCRIPTEN__
    case SDLK_v: video_toggle(); break;
#endif
    default: break;
    }
}

/* ── per-frame callback ─────────────────────────────────────────────────── */
static void loop_iter(void)
{
    SDL_Renderer *rend = app.rend;
    Sim    *sim = &g_sim;
    Camera *cam = &app.cam;

#ifdef __EMSCRIPTEN__
    poll_resize();
#endif

    /* ── events ──────────────────────────────────────────────────────── */
    SDL_Event ev;
    while(SDL_PollEvent(&ev)){
        switch(ev.type){
        case SDL_QUIT: request_quit(); break;
        case SDL_KEYDOWN: handle_key(ev.key.keysym.sym); break;
        case SDL_WINDOWEVENT:
            if(ev.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) apply_resolution();
            break;
        case SDL_MOUSEBUTTONDOWN:
            if(ev.button.button != SDL_BUTTON_LEFT) break;
            {
                double mx, my; to_logical(ev.button.x, ev.button.y, &mx, &my);
                if(app.show_help){ app.show_help = 0; break; }
                if(mx >= SIM_W){                      /* panel: pick a body */
                    int idx = (int)((my - 128) / 22);
                    if(my >= 126 && idx>=0 && idx<N_BODIES) select_body(idx);
                    break;
                }
                /* click on a body in the scene selects it */
                int best=-1; double best_d=1e9;
                for(int i=0;i<N_BODIES;i++){
                    float sx,sy; project(cam, sim->pos[i], &sx,&sy,NULL);
                    double d = hypot(mx-sx, my-sy), lim = BODIES[i].draw_r + 6;
                    if(d < lim && d < best_d){ best_d=d; best=i; }
                }
                if(best >= 0){ select_body(best); break; }
                cam->panning  = (SDL_GetModState() & KMOD_SHIFT) ? 1 : 0;
                if(cam->panning) app.follow = 0;
                cam->dragging = 1;
                cam->last_x = mx; cam->last_y = my;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            if(ev.button.button == SDL_BUTTON_LEFT){ cam->dragging = 0; cam->panning = 0; }
            break;
        case SDL_MOUSEMOTION:
            if(cam->dragging){
                double mx, my; to_logical(ev.motion.x, ev.motion.y, &mx, &my);
                double dx = mx - cam->last_x, dy = my - cam->last_y;
                cam->last_x = mx; cam->last_y = my;
                if(cam->panning){
                    cam->pan_x += dx; cam->pan_y += dy;
                } else {
                    Mat3 Ry = mat3_rot_axis((Vec3){0,1,0}, dx*0.005);
                    Mat3 Rx = mat3_rot_axis((Vec3){1,0,0}, dy*0.005);
                    cam->rot = mat3_mul(&Ry, &cam->rot);
                    cam->rot = mat3_mul(&Rx, &cam->rot);
                }
            }
            break;
        case SDL_MOUSEWHEEL:
            cam->scale *= (ev.wheel.y > 0) ? 1.12 : 0.89;
            if(cam->scale < 3.0)    cam->scale = 3.0;
            if(cam->scale > 8000.0) cam->scale = 8000.0;
            break;
        default: break;
        }
    }

    /* ── physics ─────────────────────────────────────────────────── */
    if(!app.paused){
        for(int step=0; step<app.steps_per_frame; step++){
            sim_step(sim);
            if(app.show_belt) belt_tick(sim, &g_belt);
            if(++app.trail_tick >= 3*24/(int)DT_HOURS){     /* every 3 days  */
                app.trail_tick = 0;
                for(int i=0;i<N_BODIES;i++) trail_push(sim,i);
            }
        }
        prec_sample(&app.prec, sim);
    }
    sim_energy(sim, &app.E, &app.L);
    if(app.follow) follow_center();

    /* ── render ──────────────────────────────────────────────────── */
    set_color(rend,0,0,0,255);
    SDL_RenderClear(rend);                           /* letterbox bars      */
    set_color(rend,5,5,12,255);
    SDL_Rect scene = {(int)lround(DX(0)), (int)lround(DY(0)), (int)lround(SIM_W*g_S), (int)lround(WIN_H*g_S)};
    SDL_RenderFillRect(rend, &scene);

    /* starfield (deterministic) */
    {
        unsigned s1 = 42u, s2 = 12345u;
        set_color(rend,200,200,200,80);
        for(int s=0;s<350;s++){ double x=frand(&s1)*SIM_W, y=frand(&s1)*WIN_H; dot(rend,x,y); }
        set_color(rend,255,255,255,120);
        for(int s=0;s<150;s++){ double x=frand(&s2)*SIM_W, y=frand(&s2)*WIN_H; dot(rend,x,y); }
    }

    SDL_Rect clip = scene;                           /* keep scene out of the panel */
    SDL_RenderSetClipRect(rend, &clip);

    draw_belt(rend, cam);
    for(int i=1;i<N_BODIES;i++){
        SDL_SetRenderDrawColor(rend, BODIES[i].r/3, BODIES[i].g/3, BODIES[i].b/3, 140);
        draw_orbit(rend, cam, &BODIES[i]);
    }
    for(int i=1;i<N_BODIES;i++) draw_trail(rend, cam, sim, i);
    draw_bodies(rend, cam, sim);

    SDL_RenderSetClipRect(rend, NULL);

    draw_axes(rend, cam);
    draw_text(rend, app.font_sm, "3D Solar System \xe2\x80\x93 N-body, Yoshida-4 symplectic \xc2\xb7 IC: " IC_SOURCE, 8, 6, 100,100,120);
    {
        char st[96] = "";
        if(app.paused)       strcat(st, "PAUSED  ");
        if(sim->dt < 0)      strcat(st, "\xe2\x97\x80 REVERSED  ");
        if(sim->gr)          strcat(st, "GR  ");
        if(app.follow)       strcat(st, "FOLLOW  ");
        if(app.show_belt)    strcat(st, "BELTS  ");
        if(st[0]) draw_text(rend, app.font_lg, st, 8, 24, 255,200,80);
    }

    /* FPS counter (smoothed; also feeds the panel's speed readout) */
    {
        Uint32 now=SDL_GetTicks();
        double fps_inst = 1000.0/(now - app.last_frame + 1);
        app.last_frame = now;
        app.fps = app.fps > 0 ? 0.9*app.fps + 0.1*fps_inst : fps_inst;
        double days_per_frame = app.steps_per_frame * fabs(sim->dt) * 365.25;
        char fbuf[64];
        snprintf(fbuf,64,"%.0f fps  %.2f days/frame", app.fps, days_per_frame);
        draw_text(rend,app.font_sm,fbuf, SIM_W-200, WIN_H-20, 80,80,100);
    }

#ifndef __EMSCRIPTEN__
    if(app.recording){
        if((SDL_GetTicks()/500)%2 == 0){ set_color(rend,220,30,30,255); fill_circle(rend, SIM_W-18, 14, 7); }
        draw_text(rend, app.font_sm, "REC", SIM_W-44, 6, 220,30,30);
        char rb[32]; snprintf(rb, sizeof rb, "%ds", app.rec_frames/60);
        draw_text(rend, app.font_sm, rb, SIM_W-44, 20, 180,30,30);
    }
#endif

    draw_panel(rend, sim, cam);
    if(app.show_help) draw_help(rend);

    SDL_RenderPresent(rend);

#ifndef __EMSCRIPTEN__
    maybe_snapshot();
    video_frame();
#endif
}

/* ======================================================================
   Self-tests  (solar_system --test / make test) — no window needed
   ====================================================================== */
#ifndef __EMSCRIPTEN__
static int g_tfail = 0;
static void check(int ok, const char *name, const char *detail){
    printf("  [%s] %-44s %s\n", ok ? "PASS" : "FAIL", name, detail);
    if(!ok) g_tfail++;
}

/* Fresh simulation with only the Sun and body `keep` massive (others are
   massless and simply ride along).  keep<0 keeps all masses.             */
static void test_sim(Sim *s, int keep, int gr){
    sim_init(s);
    if(keep >= 0) for(int i=1;i<N_BODIES;i++) if(i!=keep) s->mass[i] = 0.0;
    s->gr = gr;
    compute_accel(s->pos, s->vel, s->mass, s->acc, s->gr);
    sim_energy(s, &s->E0, &s->L0);
}

/* Secular rate of body i's perihelion over `years`, arcsec/century. */
static double test_precession(Sim *s, int body, double years){
    Precession p; prec_reset(&p, s, body);
    double t_end = s->time_yr + years;
    while(s->time_yr < t_end){ sim_step(s); prec_sample(&p, s); }
    double rate=0, se=0, span=0; prec_rate(&p, &rate, &se, &span);
    return rate;
}

/* Osculating semi-major axis of body i from the heliocentric energy. */
static double osc_a(const Sim *s, int i){
    Vec3 r = v3sub(s->pos[i], s->pos[0]), v = v3sub(s->vel[i], s->vel[0]);
    double mu = G_AU*(s->mass[0]+s->mass[i]);
    double eps = 0.5*v3dot(v,v) - mu/v3len(r);
    return -mu/(2.0*eps);
}

/* --dump-state: print the current initial conditions in horizons_j2000.h
   format (AU, AU/day).  Lets the Horizons-header code path be exercised
   without network access; the real header comes from tests/fetch_j2000_state.sh. */
static int dump_state(void){
    static Sim s; sim_init(&s);
    printf("/* Generated by solar_system --dump-state from %s — NOT Horizons data. */\n", IC_SOURCE);
    printf("#define HORIZONS_J2000_JD %.1f\nstatic const double HORIZONS_J2000[9][6] = {\n", JD_J2000);
    for(int i=0;i<N_BODIES;i++)
        printf("  { %.16e, %.16e, %.16e,  %.16e, %.16e, %.16e },\n",
               s.pos[i].x, s.pos[i].y, s.pos[i].z,
               s.vel[i].x/365.25, s.vel[i].y/365.25, s.vel[i].z/365.25);
    printf("};\n");
    return 0;
}

static int run_tests(void){
    static Sim s;
    char d[160];
    printf("solar_system self-tests  (dt = %g h, Yoshida-4, IC: %s)\n", DT_HOURS, IC_SOURCE);

    /* 1 ── initial state: vis-viva (table IC) / osculating a vs catalogue (Horizons IC) */
    {
        test_sim(&s, -1, 0);
        double worst = 0;
#ifdef HAVE_HORIZONS_J2000
        for(int i=1;i<N_BODIES;i++){
            double err = fabs(osc_a(&s,i) - BODIES[i].a)/BODIES[i].a; if(err > worst) worst = err;
        }
        snprintf(d, sizeof d, "max |a_osc − a_table|/a = %.1e", worst);
        check(worst < 1e-2, "Horizons vectors: frame & units sanity", d);
#else
        for(int i=1;i<N_BODIES;i++){
            Vec3 r = v3sub(s.pos[i], s.pos[0]), v = v3sub(s.vel[i], s.vel[0]);
            double mu = G_AU*(s.mass[0]+s.mass[i]);
            double v2 = v3dot(v,v), vv = mu*(2.0/v3len(r) - 1.0/BODIES[i].a);
            double err = fabs(v2-vv)/v2; if(err > worst) worst = err;
        }
        snprintf(d, sizeof d, "max |v²−μ(2/r−1/a)|/v² = %.1e", worst);
        check(worst < 1e-9, "initial state satisfies vis-viva", d);
#endif
    }

    /* 2 ── Earth two-body period = 2π√(a³/μ), a from the initial state */
    {
        test_sim(&s, 3, 0);
        double mu = G_AU*(s.mass[0]+s.mass[3]);
        double T_exp = 2*M_PI*sqrt(pow(osc_a(&s,3),3)/mu);
        Vec3 r0 = v3sub(s.pos[3], s.pos[0]);
        double th_prev = atan2(r0.y, r0.x), th_unw = 0, t_prev = 0;
        const int N = 10; double t_cross = -1;
        while(t_cross < 0){
            sim_step(&s);
            Vec3 r = v3sub(s.pos[3], s.pos[0]);
            double th = atan2(r.y, r.x), dd = th - th_prev;
            if(dd >  M_PI) dd -= 2*M_PI;
            if(dd < -M_PI) dd += 2*M_PI;
            double th_new = th_unw + dd;
            if(th_new >= N*2*M_PI){
                double f = (N*2*M_PI - th_unw)/(th_new - th_unw);
                t_cross = t_prev + f*(s.time_yr - t_prev);
            }
            th_unw = th_new; th_prev = th; t_prev = s.time_yr;
        }
        double T = t_cross/N, err = fabs(T-T_exp)/T_exp;
        snprintf(d, sizeof d, "T = %.7f yr, Kepler %.7f, rel err %.1e", T, T_exp, err);
        check(err < 1e-6, "Earth two-body period (Kepler III)", d);
    }

    /* 3 ── 9-body conservation over 100 yr */
    {
        test_sim(&s, -1, 0);
        double maxE=0, maxL=0, E, L;
        while(s.time_yr < 100.0){
            sim_step(&s);
            sim_energy(&s, &E, &L);
            double dE = fabs((E-s.E0)/s.E0), dL = fabs((L-s.L0)/s.L0);
            if(dE>maxE) maxE=dE;
            if(dL>maxL) maxL=dL;
        }
        snprintf(d, sizeof d, "max|ΔE/E₀| = %.1e, max|ΔL/L₀| = %.1e", maxE, maxL);
        check(maxE < 1e-8 && maxL < 1e-10, "energy & angular momentum, 100 yr N-body", d);
    }

    /* 4 ── Mercury: numerical drift small, GR adds 42.98″/cy */
    {
        test_sim(&s, 1, 0); double r0 = test_precession(&s, 1, 100.0);
        test_sim(&s, 1, 1); double r1 = test_precession(&s, 1, 100.0);
        snprintf(d, sizeof d, "Newtonian two-body: %+.2f″/cy", r0);
        check(fabs(r0) < 5.0, "Mercury numerical perihelion drift", d);
        snprintf(d, sizeof d, "GR − Newtonian = %.2f″/cy (theory 42.98)", r1-r0);
        check(fabs((r1-r0) - 42.98) < 1.0, "Mercury 1PN perihelion advance", d);
    }

    /* 5 ── optional: positions vs JPL Horizons (tests/horizons.csv) */
    {
        FILE *f = fopen("tests/horizons.csv", "r");
        if(!f){
            printf("  [SKIP] %-44s %s\n", "positions vs JPL Horizons",
                   "no tests/horizons.csv — run tests/fetch_horizons.sh");
        } else {
            typedef struct { int body; double jd, x, y, z; } HRow;
            HRow rows[64]; int nr = 0;
            char ln[256];
            while(fgets(ln, sizeof ln, f) && nr < 64){
                char name[32]; double jd,x,y,z;
                if(ln[0]=='#' || sscanf(ln, " %31[^,],%lf,%lf,%lf,%lf", name,&jd,&x,&y,&z) != 5) continue;
                int b=-1; for(int i=1;i<N_BODIES;i++) if(!strcasecmp(name, BODIES[i].name)) b=i;
                if(b<0) continue;
                rows[nr].body=b; rows[nr].jd=jd; rows[nr].x=x; rows[nr].y=y; rows[nr].z=z; nr++;
            }
            fclose(f);
            /* sort by JD, integrate once, landing exactly on each epoch */
            for(int a=0;a<nr;a++) for(int b=a+1;b<nr;b++) if(rows[b].jd < rows[a].jd){ HRow t=rows[a]; rows[a]=rows[b]; rows[b]=t; }
            test_sim(&s, -1, 1);          /* GR on: DE441 includes it        */
            double dt0 = s.dt; int all_ok = 1;
            for(int k=0;k<nr;k++){
                double t_target = (rows[k].jd - JD_J2000)/365.25;
                while(s.time_yr < t_target - 1e-12){
                    double rem = t_target - s.time_yr;
                    s.dt = rem < dt0 ? rem : dt0;
                    sim_step(&s);
                }
                s.dt = dt0;
                Vec3 r = v3sub(s.pos[rows[k].body], s.pos[0]);
                double err = v3len(v3sub(r, (Vec3){rows[k].x, rows[k].y, rows[k].z}));
#ifdef HAVE_HORIZONS_J2000
                /* Exact osculating start.  Achieved: ~1e-5 AU inner, ~1e-7 AU
                   outer after 20 yr; what remains is DE441 physics we omit
                   (asteroids, the Moon as a separate body, solar J2).  The
                   tolerance grows with the span so later epochs still pass. */
                double yrs = fabs(t_target); if(yrs < 20.0) yrs = 20.0;
                double tol = 1e-4*BODIES[rows[k].body].a*(yrs/20.0);
                int ok = err < tol; if(!ok) all_ok = 0;
                snprintf(d, sizeof d, "%-8s JD %.1f: |Δr| = %.2e AU (%.0f km, tol %.1e)",
                         BODIES[rows[k].body].name, rows[k].jd, err, err*1.495978707e8, tol);
                check(ok, "position vs JPL Horizons", d);
#else
                /* mean elements are not osculating: %-level errors are expected
                   (Saturn worst, via the Jupiter–Saturn great inequality), so
                   this is informational.  Generate horizons_j2000.h for a
                   strict check.                                          */
                printf("  [INFO] %-44s %-8s JD %.1f: |Δr| = %.2e AU = %.2f %% of a\n",
                       "position vs JPL Horizons (mean-element IC)",
                       BODIES[rows[k].body].name, rows[k].jd, err, 100.0*err/BODIES[rows[k].body].a);
                if(err > 0.1*BODIES[rows[k].body].a) all_ok = 0;   /* gross error only */
#endif
            }
            if(nr==0) check(0, "positions vs JPL Horizons", "file had no usable rows");
#ifndef HAVE_HORIZONS_J2000
            else check(all_ok, "positions vs JPL Horizons (gross check, < 10 % of a)",
                       "run tests/fetch_j2000_state.sh for a strict 0.1 % check");
#endif
            (void)all_ok;
        }
    }

    printf("%s\n", g_tfail ? "SOME TESTS FAILED" : "all tests passed");
    return g_tfail ? 1 : 0;
}
#endif /* !__EMSCRIPTEN__ */

/* ======================================================================
   main
   ====================================================================== */
int main(int argc, char **argv){
#ifndef __EMSCRIPTEN__
    if(argc > 1 && strcmp(argv[1], "--test") == 0)       return run_tests();
    if(argc > 1 && strcmp(argv[1], "--dump-state") == 0) return dump_state();
#else
    (void)argc; (void)argv;
#endif
    if(SDL_Init(SDL_INIT_VIDEO) < 0){ fprintf(stderr,"SDL_Init: %s\n", SDL_GetError()); return 1; }
    if(TTF_Init() < 0){ fprintf(stderr,"TTF_Init: %s\n", TTF_GetError()); return 1; }

    /* Back-buffer size.  Browser: the page measured its container in device
       pixels before the runtime started.  Native: resizable + HiDPI; the
       scene is letterboxed into whatever the drawable size is.           */
    int W = WIN_W, H = WIN_H;
    Uint32 win_flags = SDL_WINDOW_SHOWN;
#ifdef __EMSCRIPTEN__
    int dw = EM_ASM_INT({ return (window._devW | 0) || 1280; });
    if(dw < 640) dw = 640;
    if(dw > 3840) dw = 3840;
    g_last_devw = dw;
    W = dw & ~1; H = (int)lround((double)W*WIN_H/WIN_W) & ~1;
#else
    win_flags |= SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI;
    {
        const char *env = getenv("SOLAR_SCALE");     /* e.g. 1.5 → 2400×1350 */
        double sc = env ? atof(env) : 0.0;
        if(sc > 0.25 && sc < 4.0){ W = (int)lround(WIN_W*sc) & ~1; H = (int)lround(WIN_H*sc) & ~1; }
    }
#endif
    app.win = SDL_CreateWindow("Solar System", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                               W, H, win_flags);
    if(!app.win){ fprintf(stderr,"SDL_CreateWindow: %s\n", SDL_GetError()); return 1; }
    app.rend = SDL_CreateRenderer(app.win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if(!app.rend)   /* no GPU (e.g. SDL_VIDEODRIVER=dummy): software fallback */
        app.rend = SDL_CreateRenderer(app.win, -1, SDL_RENDERER_SOFTWARE);
    if(!app.rend){ fprintf(stderr,"SDL_CreateRenderer: %s\n", SDL_GetError()); return 1; }
    SDL_SetRenderDrawBlendMode(app.rend, SDL_BLENDMODE_BLEND);

    /* ── font: first path that opens wins ─────────────────────────────
       WASM: embedded in the binary's virtual FS at build time.
       Native: SOLAR_FONT, the installed copy (SOLAR_DATADIR), the source
       directory, then common system locations.                          */
#ifdef __EMSCRIPTEN__
    const char *font_paths[] = { "DejaVuSans.ttf", NULL };
#else
    const char *env_font = getenv("SOLAR_FONT");
    const char *font_paths[] = {
        env_font ? env_font : SOLAR_DATADIR "/DejaVuSans.ttf",
        SOLAR_DATADIR "/DejaVuSans.ttf",
        "DejaVuSans.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
        "/opt/homebrew/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
        NULL
    };
#endif
    for(int i=0; font_paths[i] && !app.font_lg; i++){
        TTF_Font *probe = TTF_OpenFont(font_paths[i], 12);
        if(!probe) continue;
        TTF_CloseFont(probe);
        snprintf(app.font_path, sizeof app.font_path, "%s", font_paths[i]);
        open_fonts(1.0);
    }
    if(!app.font_lg){
        fprintf(stderr, "Could not open any TTF font.\n"
#ifdef __EMSCRIPTEN__
            "Ensure DejaVuSans.ttf was embedded (--embed-file DejaVuSans.ttf).\n");
#else
            "Run `make install` (installs the bundled DejaVuSans.ttf), install DejaVu\n"
            "fonts (Debian: fonts-dejavu-core, Arch: ttf-dejavu, macOS: font-dejavu),\n"
            "or set SOLAR_FONT=/path/to/font.ttf\n");
#endif
        return 1;
    }

    /* ── init simulation & camera ────────────────────────────────────── */
    sim_init(&g_sim);
    app.cam             = cam_init();
    app.steps_per_frame = 10;                 /* 10 × 2 h = 0.83 d/frame   */
    app.running         = 1;
    app.last_frame      = SDL_GetTicks();
    app.mouse_k         = 1.0;
    prec_reset(&app.prec, &g_sim, 1);
    apply_resolution();

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(loop_iter, 0, 1);   /* never returns          */
#else
    while(app.running) loop_iter();

    if(app.recording) video_toggle();
    TTF_CloseFont(app.font_lg);
    TTF_CloseFont(app.font_sm);
    TTF_Quit();
    SDL_DestroyRenderer(app.rend);
    SDL_DestroyWindow(app.win);
    SDL_Quit();
#endif
    return 0;
}
