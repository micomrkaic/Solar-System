#define _GNU_SOURCE
/*
 * solar_system.c  –  3-D Solar System Simulation
 *
 * Physics  : Störmer-Verlet (symplectic leapfrog) N-body integrator.
 *            Conserves the symplectic structure of Hamilton's equations,
 *            giving excellent long-term energy conservation.
 *
 * Rendering: Software 3-D projection onto an SDL2 2-D renderer.
 *            Arcball (trackball) rotation via mouse drag.
 *            Zoom via scroll wheel.
 *
 * Units    : AU (astronomical units), solar masses, years.
 *            G = 4π²  in these units.
 *
 * Compile  : gcc -O2 -std=c17 solar_system.c -o solar_system \
 *                $(pkg-config --cflags --libs sdl2 SDL2_ttf)
 *
 * Recording: Press V to start/stop. Requires ffmpeg in PATH.
 *            Output: solar_YYYYMMDD_HHMMSS.mp4
 *
 * SDL3 port: Replace SDL_RenderDrawPointF → SDL_RenderPoint,
 *            SDL_RenderDrawLineF → SDL_RenderLine, and update
 *            the init/destroy boilerplate.  Everything else is identical.
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ── constants ─────────────────────────────────────────────────────────── */
#define G_AU        (4.0 * M_PI * M_PI)   /* AU³/(M☉·yr²)               */
#define N_BODIES    9                       /* Sun + 8 planets             */
#define TRAIL_LEN   21000                   /* trail history points (fits Neptune) */
#define WIN_W       1600
#define WIN_H       900
#define PANEL_W     310                     /* right-side info panel width */
#define SIM_W       (WIN_W - PANEL_W)

/* ── body data ─────────────────────────────────────────────────────────── */
typedef struct {
    const char *name;
    double mass;            /* solar masses                                 */
    double a;               /* semi-major axis [AU]                         */
    double ecc;             /* eccentricity                                 */
    double inc;             /* inclination [deg]                            */
    double omega;           /* longitude of ascending node [deg]            */
    double w;               /* argument of perihelion [deg]                 */
    double M0;              /* mean anomaly at J2000 [deg]                  */
    double radius_km;       /* physical radius [km]  (display only)         */
    double period_yr;       /* sidereal period [yr]  (display only)         */
    double v_orb_km_s;      /* mean orbital speed [km/s] (display only)     */
    Uint8  r, g, b;         /* display colour                               */
    double draw_r;          /* on-screen radius [pixels]                    */
} BodyDef;

/* Orbital elements from JPL Horizons (J2000.0 epoch).
   Masses: Sun=1, others in solar masses. */
static const BodyDef BODIES[N_BODIES] = {
/* All orbital elements: JPL SSD Table 1, J2000.0 epoch
   ω = long.peri − long.node,  M0 = mean.lon − long.peri  (both mod 360°)
   name        mass(M☉)      a(AU)       ecc        inc(°)     Ω(°)        ω(°)       M0(°)      R(km)     T(yr)  v(km/s) R   G   B   dr */
  { "Sun",    1.0,          0.0,        0.0,        0.0,       0.0,        0.0,        0.0,   695700.0,    0.0,    0.0, 255,220, 50, 14.0 },
  { "Mercury",1.6601e-7,    0.38709927, 0.20563593, 7.00498,  48.33077,   29.12703,  174.79,   2439.7,   0.2408,  47.4, 180,180,180,  3.0 },
  { "Venus",  2.4478e-6,    0.72333566, 0.00677672, 3.39468,  76.67984,   54.92262,   50.38,   6051.8,   0.6152,  35.0, 230,180, 80,  5.0 },
  { "Earth",  3.0034e-6,    1.00000261, 0.01671123, 0.00000,   0.0,       102.93768,  357.53,  6371.0,   1.0000,  29.8,  80,130,255,  5.0 },
  { "Mars",   3.2271e-7,    1.52371034, 0.09339410, 1.84969,  49.55954,  286.49683,   19.39,   3389.5,   1.8809,  24.1, 200, 80, 50,  4.0 },
  { "Jupiter",9.5458e-4,    5.20288700, 0.04838624, 1.30440, 100.47391,  274.25457,   19.65,  69911.0,  11.8618,  13.1, 220,160, 90, 11.0 },
  { "Saturn", 2.8580e-4,    9.53667594, 0.05386179, 2.48599, 113.66242,  338.93646,  317.36,  58232.0,  29.4571,   9.7, 210,190,130, 10.0 },
  { "Uranus", 4.3653e-5,   19.18916464, 0.04725744, 0.77264,  74.01693,   96.93735,  142.28,  25362.0,  84.0107,   6.8, 150,220,220,  7.0 },
  { "Neptune",5.1514e-5,   30.06992276, 0.00859048, 1.77004, 131.78423,  273.18054,  259.92,  24622.0, 164.7945,   5.4,  80,100,255,  7.0 },
};

/* ── vector helpers ─────────────────────────────────────────────────────── */
typedef struct { double x, y, z; } Vec3;

static inline Vec3 v3add(Vec3 a, Vec3 b) { return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 v3sub(Vec3 a, Vec3 b) { return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 v3scale(Vec3 a, double s){ return (Vec3){a.x*s, a.y*s, a.z*s}; }
static inline double v3dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double v3len(Vec3 a) { return sqrt(v3dot(a,a)); }
static inline Vec3 v3norm(Vec3 a)  { double l=v3len(a); return l>0?v3scale(a,1.0/l):(Vec3){0,0,0}; }
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

/* ── simulation state ───────────────────────────────────────────────────── */
typedef struct {
    Vec3 pos[N_BODIES];
    Vec3 vel[N_BODIES];       /* leapfrog half-step velocities              */
    double mass[N_BODIES];

    /* trail ring buffer */
    Vec3  trail[N_BODIES][TRAIL_LEN];
    int   trail_head[N_BODIES];
    int   trail_count[N_BODIES];
    int   trail_cap[N_BODIES];  /* per-body capacity (≤ TRAIL_LEN)        */

    double time_yr;           /* elapsed simulation time [yr]               */
    double dt;                /* time step [yr]                             */
    double speed;             /* simulation speed multiplier                */
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

/* convert orbital elements → Cartesian state (in ecliptic J2000) */
static void elements_to_cartesian(const BodyDef *bd, Vec3 *pos, Vec3 *vel){
    if(bd->a == 0.0){ *pos=(Vec3){0,0,0}; *vel=(Vec3){0,0,0}; return; }

    double a   = bd->a;
    double e   = bd->ecc;
    double inc = bd->inc   * M_PI/180.0;
    double Om  = bd->omega * M_PI/180.0;
    double w   = bd->w    * M_PI/180.0;
    double M0  = bd->M0   * M_PI/180.0;

    double E   = solve_kepler(M0, e);

    /* position & velocity in orbital plane */
    double nu  = 2.0*atan2(sqrt(1+e)*sin(E/2.0), sqrt(1-e)*cos(E/2.0));
    double r   = a*(1 - e*cos(E));
    double px  = r*cos(nu), py = r*sin(nu);

    double n   = 2.0*M_PI / bd->period_yr;           /* mean motion [rad/yr] */
    /* vx, vy in orbital frame — correct Keplerian formula (Bate, Mueller & White):
       vfac = n*a^2/r  (NOT a*n/sqrt(1-e^2), which is only correct at semi-latus rectum) */
    double vfac = n*a*a / r;
    double vx  = -vfac*sin(E);
    double vy  =  vfac*sqrt(1-e*e)*cos(E);

    /* rotate to ecliptic: R3(-Om) · R1(-inc) · R3(-w) */
    double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);

    /* Perifocal → ecliptic rotation matrix elements */
    double Qxx = cO*cw - sO*sw*ci,  Qxy = -cO*sw - sO*cw*ci;
    double Qyx = sO*cw + cO*sw*ci,  Qyy = -sO*sw + cO*cw*ci;
    double Qzx = sw*si,              Qzy =  cw*si;

    pos->x = Qxx*px + Qxy*py;
    pos->y = Qyx*px + Qyy*py;
    pos->z = Qzx*px + Qzy*py;

    vel->x = Qxx*vx + Qxy*vy;
    vel->y = Qyx*vx + Qyy*vy;
    vel->z = Qzx*vx + Qzy*vy;
}

/* ── compute gravitational accelerations ────────────────────────────────── */
static void compute_accel(const Vec3 *pos, const double *mass, Vec3 *acc){
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
}

/* ── Störmer-Verlet leapfrog step ────────────────────────────────────────── */
static void leapfrog_step(Sim *sim){
    Vec3 acc[N_BODIES];
    compute_accel(sim->pos, sim->mass, acc);

    double hdt = 0.5 * sim->dt;

    /* kick: v += a * dt/2 */
    for(int i=0;i<N_BODIES;i++)
        sim->vel[i] = v3add(sim->vel[i], v3scale(acc[i], hdt));

    /* drift: x += v * dt */
    for(int i=0;i<N_BODIES;i++)
        sim->pos[i] = v3add(sim->pos[i], v3scale(sim->vel[i], sim->dt));

    /* recompute accel at new positions */
    compute_accel(sim->pos, sim->mass, acc);

    /* kick: v += a * dt/2 */
    for(int i=0;i<N_BODIES;i++)
        sim->vel[i] = v3add(sim->vel[i], v3scale(acc[i], hdt));

    sim->time_yr += sim->dt;
}

/* ── simulation init ─────────────────────────────────────────────────────── */
static void sim_init(Sim *sim){
    memset(sim,0,sizeof(*sim));
    sim->dt    = 1.0/(365.25*6.0); /* 1/6 day in years — finer minimum speed  */
    sim->speed = 1.0;

    /* First pass: set all masses, positions, velocities from Keplerian elements */
    for(int i=0;i<N_BODIES;i++){
        sim->mass[i] = BODIES[i].mass;
        elements_to_cartesian(&BODIES[i], &sim->pos[i], &sim->vel[i]);
    }
    /* Second pass: single correct leapfrog half-kick using full N-body accel */
    {
        Vec3 acc[N_BODIES];
        compute_accel(sim->pos, sim->mass, acc);
        for(int i=0;i<N_BODIES;i++)
            sim->vel[i] = v3sub(sim->vel[i], v3scale(acc[i], 0.5*sim->dt));
    }
    /* compute per-body trail capacities: one full orbit each */
    for(int i=0;i<N_BODIES;i++){
        if(BODIES[i].period_yr <= 0.0){
            sim->trail_cap[i] = 0;   /* Sun: no trail */
        } else {
            /* 1 sample every 3 physics steps, each step = 1 day */
            double samples_per_yr = 365.25 / 3.0;
            int cap = (int)(BODIES[i].period_yr * samples_per_yr);
            if(cap < 50)       cap = 50;
            if(cap > TRAIL_LEN) cap = TRAIL_LEN;
            sim->trail_cap[i] = cap;
        }
    }
}

/* ── trail management ───────────────────────────────────────────────────── */
static void trail_push(Sim *sim, int i){
    int cap = sim->trail_cap[i];
    if(cap <= 0) return;           /* Sun: skip */
    int h = sim->trail_head[i];
    sim->trail[i][h] = sim->pos[i];
    sim->trail_head[i] = (h+1) % cap;
    if(sim->trail_count[i] < cap) sim->trail_count[i]++;
}

/* ── camera / view state ─────────────────────────────────────────────────── */
typedef struct {
    Mat3   rot;          /* current rotation matrix (arcball)               */
    double scale;        /* pixels per AU                                   */
    double fov_dist;     /* "focal length" for perspective                  */
    int    dragging;
    int    panning;      /* shift+drag pan mode                             */
    int    last_x, last_y;
    double pan_x, pan_y; /* screen-space translation offset [pixels]       */
    int    selected;     /* which body is selected (-1=none)                */
} Camera;

static Camera cam_init(void){
    Camera c;
    c.rot     = mat3_identity();
    c.scale   = 60.0;     /* inner system fit                               */
    c.fov_dist = 8.0;     /* perspective foreshortening                     */
    c.dragging = 0;
    c.panning  = 0;
    c.last_x   = 0;
    c.last_y   = 0;
    c.pan_x    = 0.0;
    c.pan_y    = 0.0;
    c.selected = 3;       /* start with Earth selected                      */
    /* tilt slightly to show 3-D structure */
    Vec3 xax = {1,0,0};
    Mat3 tilt = mat3_rot_axis(xax, -0.4);
    c.rot = mat3_mul(&tilt, &c.rot);
    return c;
}


/* Set camera to a canonical axis-aligned view, preserving pan/zoom/selection */
static void cam_set_view(Camera *c, int axis){
    /* axis: 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane / top-down) */
    Mat3 m = mat3_identity();
    Vec3 ax;
    double angle;
    switch(axis){
        case 0:  /* look down +X → rotate Y-axis -90° */
            ax = (Vec3){0,1,0}; angle = -M_PI/2.0;
            m = mat3_rot_axis(ax, angle);
            break;
        case 1:  /* look down +Y → rotate X-axis +90° */
            ax = (Vec3){1,0,0}; angle = M_PI/2.0;
            m = mat3_rot_axis(ax, angle);
            break;
        case 2:  /* look down +Z → identity (XY plane, top-down) */
            m = mat3_identity();
            break;
    }
    c->rot = m;
}

/* project a 3-D point to screen */
static void project(const Camera *cam, Vec3 world,
                    float *sx, float *sy, float *sz_out){
    Vec3 v = mat3_mul_vec(&cam->rot, world);
    /* perspective divide */
    double z_ofs = cam->fov_dist;
    double w     = z_ofs / (z_ofs + v.z + 40.0);   /* +40 keeps w > 0    */
    *sx = (float)(SIM_W/2.0 + cam->pan_x + v.x * cam->scale * w);
    *sy = (float)(WIN_H/2.0 + cam->pan_y - v.y * cam->scale * w);
    if(sz_out) *sz_out = (float)v.z;
}

/* ── SDL helpers ─────────────────────────────────────────────────────────── */
static void set_color(SDL_Renderer *rend, Uint8 r, Uint8 g, Uint8 b, Uint8 a){
    SDL_SetRenderDrawColor(rend, r, g, b, a);
}

/* draw a filled circle (midpoint algorithm) */
static void fill_circle(SDL_Renderer *rend, int cx, int cy, int rad){
    for(int dy=-rad; dy<=rad; dy++){
        int dx = (int)sqrt((double)(rad*rad - dy*dy));
        SDL_RenderDrawLine(rend, cx-dx, cy+dy, cx+dx, cy+dy);
    }
}

/* draw text (simple wrapper) */
static void draw_text(SDL_Renderer *rend, TTF_Font *font,
                      const char *text, int x, int y,
                      Uint8 r, Uint8 g, Uint8 b){
    SDL_Color col = {r,g,b,255};
    SDL_Surface *surf = TTF_RenderUTF8_Blended(font, text, col);
    if(!surf) return;
    SDL_Texture *tex = SDL_CreateTextureFromSurface(rend, surf);
    if(tex){
        SDL_Rect dst = {x, y, surf->w, surf->h};
        SDL_RenderCopy(rend, tex, NULL, &dst);
        SDL_DestroyTexture(tex);
    }
    SDL_FreeSurface(surf);
}


/* ── render orbit ellipse ───────────────────────────────────────────────── */
static void draw_orbit(SDL_Renderer *rend, const Camera *cam,
                       const BodyDef *bd, Uint8 r, Uint8 g, Uint8 b){
    if(bd->a == 0.0) return;
    int STEPS = 180;
    float px=0, py=0;
    for(int step=0; step<=STEPS; step++){
        double nu = 2.0*M_PI*step/STEPS;
        double dist = bd->a*(1.0 - bd->ecc*bd->ecc)/(1.0 + bd->ecc*cos(nu));

        /* position in orbital plane */
        double ox = dist*cos(nu), oy = dist*sin(nu);

        /* rotate to ecliptic (same as elements_to_cartesian) */
        double inc = bd->inc   * M_PI/180.0;
        double Om  = bd->omega * M_PI/180.0;
        double w   = bd->w    * M_PI/180.0;
        double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);
        double Qxx = cO*cw - sO*sw*ci,  Qxy = -cO*sw - sO*cw*ci;
        double Qyx = sO*cw + cO*sw*ci,  Qyy = -sO*sw + cO*cw*ci;
        double Qzx = sw*si,              Qzy =  cw*si;
        Vec3 p = {Qxx*ox+Qxy*oy, Qyx*ox+Qyy*oy, Qzx*ox+Qzy*oy};

        float sx,sy; project(cam, p, &sx, &sy, NULL);
        if(step>0)
            SDL_RenderDrawLineF(rend, px, py, sx, sy);
        px=sx; py=sy;
    }
}

/* ── draw the info panel ───────────────────────────────────────────────── */
static void draw_panel(SDL_Renderer *rend, TTF_Font *font_lg, TTF_Font *font_sm,
                       const Sim *sim, const Camera *cam, int sel){
    int px0 = SIM_W + 8;

    /* panel background */
    set_color(rend,15,15,25,255);
    SDL_Rect bg = {SIM_W, 0, PANEL_W, WIN_H};
    SDL_RenderFillRect(rend, &bg);

    /* separator line */
    set_color(rend,60,60,80,255);
    SDL_RenderDrawLine(rend, SIM_W, 0, SIM_W, WIN_H);

    /* title */
    draw_text(rend, font_lg, "Solar System", px0, 8, 255,220,50);

    char buf[128];
    /* elapsed time */
    double yr = sim->time_yr;
    int years = (int)yr;
    int days  = (int)((yr - years)*365.25);
    snprintf(buf, sizeof(buf), "T = %d yr %d d", years, days);
    draw_text(rend, font_sm, buf, px0, 36, 180,180,200);

    /* speed */
    snprintf(buf, sizeof(buf), "Speed: %.0f× real", sim->speed * sim->dt / (1.0/365.25));
    draw_text(rend, font_sm, buf, px0, 52, 180,180,200);

    /* divider */
    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, 72, WIN_W-4, 72);

    /* planet list  – highlight selected */
    int y0 = 80;
    for(int i=0;i<N_BODIES;i++){
        const BodyDef *bd = &BODIES[i];
        if(i == sel){
            SDL_Rect hl = {SIM_W, y0-2, PANEL_W, 22};
            set_color(rend,30,30,55,255);
            SDL_RenderFillRect(rend, &hl);
        }
        set_color(rend, bd->r, bd->g, bd->b, 255);
        fill_circle(rend, px0+6, y0+8, 5);
        draw_text(rend, font_sm, bd->name, px0+16, y0,
                  (i==sel)?255:200, (i==sel)?255:200, (i==sel)?255:200);
        y0 += 22;
    }

    /* divider */
    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, y0+4, WIN_W-4, y0+4);

    /* selected body details */
    if(sel >= 0 && sel < N_BODIES){
        const BodyDef *bd = &BODIES[sel];
        int yd = y0 + 14;
        draw_text(rend, font_lg, bd->name, px0, yd, bd->r, bd->g, bd->b);
        yd += 26;

        /* distance from Sun */
        Vec3 r = v3sub(sim->pos[sel], sim->pos[0]);
        double dist_AU = v3len(r);
        double dist_km = dist_AU * 1.496e8;

        char vals[20][64];
        vals[0][0] = '\0';              /* section header */
        vals[1][0] = '\0';
        snprintf(vals[2], 64, "%.4f AU",  bd->a);
        vals[3][0] = '\0';
        snprintf(vals[4], 64, "%.5f",     bd->ecc);
        vals[5][0] = '\0';
        snprintf(vals[6], 64, "%.3f°",    bd->inc);
        vals[7][0] = '\0';
        if(bd->period_yr > 0)
            snprintf(vals[8], 64, "%.3f yr", bd->period_yr);
        else
            snprintf(vals[8], 64, "–");
        vals[9][0] = '\0';
        vals[10][0] = '\0';             /* spacer */
        vals[11][0] = '\0';
        vals[12][0] = '\0';             /* section header */
        vals[13][0] = '\0';
        snprintf(vals[14], 64, "%.4f AU (%.2e km)", dist_AU, dist_km);
        vals[15][0] = '\0';

        /* current speed magnitude */
        Vec3 vrel = v3sub(sim->vel[sel], sim->vel[0]);
        double v_AU_yr = v3len(vrel);
        /* 1 AU/yr = 4.740 km/s */
        double v_km_s = v_AU_yr * 4.740;
        snprintf(vals[16], 64, "%.2f km/s", v_km_s);
        vals[17][0] = '\0';
        snprintf(vals[18], 64, "%.2f km/s", bd->v_orb_km_s);
        vals[19][0] = '\0';

        /* Row titles and values */
        struct { const char *lbl; const char *val; int section; } rows[] = {
            {"Orbital elements", "",    1},
            {"  Semi-major axis", vals[2],  0},
            {"  Eccentricity",   vals[4],  0},
            {"  Inclination",    vals[6],  0},
            {"  Period",         vals[8],  0},
            {"Current state",    "",       1},
            {"  Dist. from Sun", "",       0},
            {"  ",               vals[14], 0},
            {"  Speed now",      vals[16], 0},
            {"  Mean speed",     vals[18], 0},
            {"Physical",         "",       1},
            {"  Radius",         NULL,     0},
            {"  Mass (M\xe2\x8a\x99)", NULL, 0},
        };

        char rbuf[64], mbuf[64];
        snprintf(rbuf, 64, "%.0f km", bd->radius_km);
        snprintf(mbuf, 64, "%.3e M\xe2\x8a\x99", bd->mass);
        rows[11].val = rbuf;
        rows[12].val = mbuf;

        for(int k=0; k<(int)(sizeof(rows)/sizeof(rows[0])); k++){
            if(rows[k].section){
                draw_text(rend, font_sm, rows[k].lbl, px0, yd,
                          180,180,100);
            } else {
                draw_text(rend, font_sm, rows[k].lbl, px0, yd,
                          140,140,160);
                if(rows[k].val && rows[k].val[0])
                    draw_text(rend, font_sm, rows[k].val, px0+8, yd+13,
                              220,220,240);
            }
            yd += (rows[k].val && rows[k].val[0] && !rows[k].section)? 28 : 18;
        }
    }

    /* controls hint at bottom */
    int yh = WIN_H - 100;
    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, yh-4, WIN_W-4, yh-4);
    draw_text(rend, font_sm, "Controls",             px0, yh,    160,160,100);
    draw_text(rend, font_sm, "Drag        – rotate",  px0, yh+14, 120,120,140);
    draw_text(rend, font_sm, "Shift+drag  – pan",     px0, yh+28, 120,120,140);
    draw_text(rend, font_sm, "Scroll      – zoom",    px0, yh+42, 120,120,140);
    draw_text(rend, font_sm, "+/-  – faster/slower",  px0, yh+56, 120,120,140);
    draw_text(rend, font_sm, "1-9  – select body",    px0, yh+70, 120,120,140);
    draw_text(rend, font_sm, "R    – reset view",     px0, yh+84, 120,120,140);
    draw_text(rend, font_sm, "X/Y/Z – axis views",    px0, yh+98, 120,120,140);
    draw_text(rend, font_sm, "V    – record video",   px0, yh+112, 120,120,140);
    draw_text(rend, font_sm, "Q / Esc  – quit",       px0, yh+126, 120,120,140);
}

/* ── draw legend at top-left of sim area ───────────────────────────────── */
static void draw_legend(SDL_Renderer *rend, TTF_Font *font,
                        double sim_yr, double speed_mult){
    (void)sim_yr; (void)speed_mult;
    draw_text(rend, font, "3D Solar System — Störmer-Verlet", 8, 6, 100,100,120);
}

/* ── main ───────────────────────────────────────────────────────────────── */
int main(void){
    /* init SDL */
    if(SDL_Init(SDL_INIT_VIDEO) < 0){
        fprintf(stderr,"SDL_Init: %s\n", SDL_GetError()); return 1;
    }
    if(TTF_Init() < 0){
        fprintf(stderr,"TTF_Init: %s\n", TTF_GetError()); return 1;
    }

    SDL_Window   *win  = SDL_CreateWindow("Solar System",
                             SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                             WIN_W, WIN_H, SDL_WINDOW_SHOWN|SDL_WINDOW_RESIZABLE);
    SDL_Renderer *rend = SDL_CreateRenderer(win, -1,
                             SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    if(!win||!rend){ fprintf(stderr,"Window/Renderer: %s\n",SDL_GetError()); return 1; }

    /* try to load a system font; fall back to a built-in path */
    const char *font_paths[] = {
        /* macOS system fonts */
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
        /* macOS Homebrew */
        "/opt/homebrew/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/dejavu/DejaVuSans.ttf",
        /* Linux Debian/Ubuntu */
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        /* Linux Arch/Fedora */
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
        NULL
    };
    TTF_Font *font_lg = NULL, *font_sm = NULL;
    for(int i=0; font_paths[i]; i++){
        font_lg = TTF_OpenFont(font_paths[i], 15);
        font_sm = TTF_OpenFont(font_paths[i], 12);
        if(font_lg && font_sm) break;
        if(font_lg){ TTF_CloseFont(font_lg); font_lg=NULL; }
        if(font_sm){ TTF_CloseFont(font_sm); font_sm=NULL; }
    }
    if(!font_lg || !font_sm){
        fprintf(stderr,
            "Could not open any TTF font.\n"
            "Set SOLAR_FONT=/path/to/font.ttf and rerun, or install DejaVu fonts:\n"
            "  macOS : brew install font-dejavu           (needs homebrew-cask-fonts)\n"
            "  Debian: sudo apt install fonts-dejavu-core\n"
            "  Arch  : sudo pacman -S ttf-dejavu\n");
        const char *env_font = getenv("SOLAR_FONT");
        if(env_font){
            font_lg = TTF_OpenFont(env_font, 15);
            font_sm = TTF_OpenFont(env_font, 12);
        }
        if(!font_lg || !font_sm) return 1;
    }

    static Sim sim;
    sim_init(&sim);
    Camera cam = cam_init();

    /* trail sampling counter */
    int trail_tick = 0;

    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);

    int running = 1;
    SDL_Event ev;
    Uint32 last_frame = SDL_GetTicks();

    /* number of physics steps per rendered frame */
    int steps_per_frame = 5;

    /* ── video recording state ───────────────────────────────────── */
    FILE  *ffmpeg_pipe  = NULL;   /* popen'd ffmpeg process           */
    Uint8 *pixel_buf    = NULL;   /* frame pixel buffer               */
    int    recording    = 0;      /* currently recording?             */
    int    rec_frames   = 0;      /* frames written this recording    */


    while(running){
        /* ── events ──────────────────────────────────────────────────── */
        while(SDL_PollEvent(&ev)){
            if(ev.type == SDL_QUIT) running=0;
            if(ev.type == SDL_KEYDOWN){
                switch(ev.key.keysym.sym){
                case SDLK_ESCAPE: running=0; break;
                case SDLK_q: running=0; break;
                case SDLK_EQUALS: case SDLK_PLUS:
                    steps_per_frame = steps_per_frame + 1;
                    if(steps_per_frame>3000) steps_per_frame=3000;
                    break;
                case SDLK_MINUS:
                    steps_per_frame = steps_per_frame - 1;
                    if(steps_per_frame<1) steps_per_frame=1;
                    break;
                case SDLK_r: cam = cam_init(); break;
                case SDLK_v: {
                    if(!recording){
                        /* build filename with timestamp */
                        time_t t = time(NULL);
                        struct tm *tm = localtime(&t);
                        char fname[64];
                        strftime(fname, sizeof(fname),
                                 "solar_%Y%m%d_%H%M%S.mp4", tm);
                        char cmd[256];
                        snprintf(cmd, sizeof(cmd),
                            "ffmpeg -y -f rawvideo -pixel_format bgra "
                            "-video_size %dx%d -framerate 60 -i pipe:0 "
                            "-c:v libx264 -preset fast -crf 18 "
                            "-pix_fmt yuv420p %s",
                            WIN_W, WIN_H, fname);
                        ffmpeg_pipe = popen(cmd, "w");
                        if(ffmpeg_pipe){
                            pixel_buf = malloc(WIN_W * WIN_H * 4);
                            recording = 1; rec_frames = 0;
                            SDL_Log("Recording started: %s", fname);
                        } else {
                            SDL_Log("Failed to open ffmpeg pipe");
                        }
                    } else {
                        /* stop recording */
                        recording = 0;
                        free(pixel_buf); pixel_buf = NULL;
                        pclose(ffmpeg_pipe); ffmpeg_pipe = NULL;
                        SDL_Log("Recording stopped (%d frames)", rec_frames);
                    }
                } break;
                case SDLK_x: cam_set_view(&cam, 0); break; /* YZ plane */
                case SDLK_y: cam_set_view(&cam, 1); break; /* XZ plane */
                case SDLK_z: cam_set_view(&cam, 2); break; /* XY plane – top-down */
                case SDLK_1: cam.selected=0; break;
                case SDLK_2: cam.selected=1; break;
                case SDLK_3: cam.selected=2; break;
                case SDLK_4: cam.selected=3; break;
                case SDLK_5: cam.selected=4; break;
                case SDLK_6: cam.selected=5; break;
                case SDLK_7: cam.selected=6; break;
                case SDLK_8: cam.selected=7; break;
                case SDLK_9: cam.selected=8; break;
                }
            }
            if(ev.type == SDL_MOUSEBUTTONDOWN && ev.button.button==SDL_BUTTON_LEFT){
                SDL_Keymod mod = SDL_GetModState();
                cam.panning  = (mod & KMOD_SHIFT) ? 1 : 0;
                cam.dragging = 1;
                cam.last_x   = ev.button.x;
                cam.last_y   = ev.button.y;
            }
            if(ev.type == SDL_MOUSEBUTTONUP && ev.button.button==SDL_BUTTON_LEFT){
                cam.dragging = 0;
                cam.panning  = 0;
            }
            if(ev.type == SDL_MOUSEMOTION && cam.dragging){
                int dx = ev.motion.x - cam.last_x;
                int dy = ev.motion.y - cam.last_y;
                cam.last_x = ev.motion.x; cam.last_y = ev.motion.y;
                if(cam.panning){
                    /* shift+drag: translate (pan) */
                    cam.pan_x += dx;
                    cam.pan_y += dy;
                } else {
                    /* plain drag: arcball rotate around world Y then X */
                    Vec3 yax={0,1,0}, xax={1,0,0};
                    Mat3 Ry = mat3_rot_axis(yax, dx*0.005);
                    Mat3 Rx = mat3_rot_axis(xax, dy*0.005);
                    cam.rot = mat3_mul(&Ry, &cam.rot);
                    cam.rot = mat3_mul(&Rx, &cam.rot);
                }
            }
            if(ev.type == SDL_MOUSEWHEEL){
                cam.scale *= (ev.wheel.y > 0) ? 1.12 : 0.89;
                if(cam.scale<3.0)   cam.scale=3.0;
                if(cam.scale>8000.0) cam.scale=8000.0;
            }
            /* click in panel selects body */
            if(ev.type==SDL_MOUSEBUTTONDOWN && ev.button.button==SDL_BUTTON_LEFT){
                int mx=ev.button.x, my=ev.button.y;
                if(mx >= SIM_W){
                    int idx = (my - 80) / 22;
                    if(idx>=0 && idx<N_BODIES) cam.selected=idx;
                }
            }
        }

        /* ── physics ─────────────────────────────────────────────────── */
        for(int step=0; step<steps_per_frame; step++){
            leapfrog_step(&sim);
            trail_tick++;
            if(trail_tick >= 18){    /* record trail every 18 steps = 3 days */
                trail_tick=0;
                for(int i=0;i<N_BODIES;i++) trail_push(&sim,i);
            }
        }

        /* ── render ──────────────────────────────────────────────────── */
        set_color(rend,5,5,12,255);
        SDL_RenderClear(rend);

        /* starfield (static) */
        srand(42);
        set_color(rend,200,200,200,80);
        for(int s=0;s<350;s++){
            int sx=rand()%SIM_W, sy=rand()%WIN_H;
            SDL_RenderDrawPoint(rend,sx,sy);
        }
        srand(12345);
        set_color(rend,255,255,255,120);
        for(int s=0;s<150;s++){
            int sx=rand()%SIM_W, sy=rand()%WIN_H;
            SDL_RenderDrawPoint(rend,sx,sy);
        }

        /* orbit ellipses (faint) */
        for(int i=1;i<N_BODIES;i++){
            Uint8 r=BODIES[i].r, g=BODIES[i].g, b=BODIES[i].b;
            SDL_SetRenderDrawColor(rend, r/3, g/3, b/3, 80);
            draw_orbit(rend, &cam, &BODIES[i], r/3, g/3, b/3);
        }

        /* trails */
        for(int i=1;i<N_BODIES;i++){
            int cnt = sim.trail_count[i];
            int head = sim.trail_head[i];
            int cap  = sim.trail_cap[i];
            if(cap <= 0) continue;
            float px2=0,py2=0;
            for(int k=0;k<cnt;k++){
                int idx = ((head - cnt + k) + cap) % cap;
                float sx,sy;
                project(&cam, sim.trail[i][idx], &sx, &sy, NULL);
                float alpha = 180.0f * k / cnt;
                SDL_SetRenderDrawColor(rend,
                    BODIES[i].r, BODIES[i].g, BODIES[i].b, (Uint8)alpha);
                if(k>0) SDL_RenderDrawLineF(rend,px2,py2,sx,sy);
                px2=sx; py2=sy;
            }
        }

        /* bodies */
        for(int i=0;i<N_BODIES;i++){
            float sx,sy;
            project(&cam, sim.pos[i], &sx, &sy, NULL);
            int rad = (int)BODIES[i].draw_r;

            /* glow for Sun */
            if(i==0){
                for(int glow=rad+8; glow>rad; glow--){
                    Uint8 a=(Uint8)(30*(glow-rad));
                    SDL_SetRenderDrawColor(rend,255,200,50,a);
                    /* draw circle approximation */
                    for(int dy2=-glow;dy2<=glow;dy2++){
                        int dx2=(int)sqrt((double)(glow*glow-dy2*dy2));
                        SDL_RenderDrawLine(rend,(int)sx-dx2,(int)sy+dy2,(int)sx+dx2,(int)sy+dy2);
                    }
                }
            }

            set_color(rend, BODIES[i].r, BODIES[i].g, BODIES[i].b, 255);
            fill_circle(rend, (int)sx, (int)sy, rad);

            /* selection ring */
            if(i == cam.selected){
                SDL_SetRenderDrawColor(rend,255,255,255,160);
                int rr=rad+4;
                for(int ang=0;ang<360;ang+=4){
                    float a1=ang*M_PI/180.0f, a2=(ang+4)*M_PI/180.0f;
                    SDL_RenderDrawLineF(rend,
                        sx+rr*cosf(a1), sy+rr*sinf(a1),
                        sx+rr*cosf(a2), sy+rr*sinf(a2));
                }
            }

            /* label */
            draw_text(rend, font_sm, BODIES[i].name,
                      (int)sx+rad+4, (int)sy-7,
                      BODIES[i].r, BODIES[i].g, BODIES[i].b);
        }

        /* axis indicator (bottom-left of sim area) */
        {
            int ax0=50, ay0=WIN_H-50;
            float len=30.0f;
            Vec3 ox={1,0,0},oy={0,1,0},oz={0,0,1};
            Vec3 rx=mat3_mul_vec(&cam.rot,ox);
            Vec3 ry=mat3_mul_vec(&cam.rot,oy);
            Vec3 rz=mat3_mul_vec(&cam.rot,oz);
            SDL_SetRenderDrawColor(rend,200,50,50,200);
            SDL_RenderDrawLineF(rend,ax0,ay0, ax0+len*rx.x, ay0-len*rx.y);
            draw_text(rend,font_sm,"X",ax0+(int)(len*rx.x)+2,ay0-(int)(len*rx.y),200,50,50);
            SDL_SetRenderDrawColor(rend,50,200,50,200);
            SDL_RenderDrawLineF(rend,ax0,ay0, ax0+len*ry.x, ay0-len*ry.y);
            draw_text(rend,font_sm,"Y",ax0+(int)(len*ry.x)+2,ay0-(int)(len*ry.y),50,200,50);
            SDL_SetRenderDrawColor(rend,50,100,220,200);
            SDL_RenderDrawLineF(rend,ax0,ay0, ax0+len*rz.x, ay0-len*rz.y);
            draw_text(rend,font_sm,"Z",ax0+(int)(len*rz.x)+2,ay0-(int)(len*rz.y),50,100,220);
        }

        draw_legend(rend, font_sm, sim.time_yr, (double)steps_per_frame);

        /* FPS counter */
        {
            Uint32 now=SDL_GetTicks();
            double fps = 1000.0/(now-last_frame+1);
            last_frame=now;
            double days_per_frame = steps_per_frame * sim.dt * 365.25;
            char fbuf[64]; snprintf(fbuf,64,"%.0f fps  %.2f days/frame", fps, days_per_frame);
            draw_text(rend,font_sm,fbuf, SIM_W-200, WIN_H-20, 80,80,100);
        }

        /* REC indicator */
        if(recording){
            /* blinking red dot */
            Uint32 ms = SDL_GetTicks();
            if((ms/500)%2 == 0){  /* blink every 500ms */
                set_color(rend, 220, 30, 30, 255);
                fill_circle(rend, SIM_W - 18, 14, 7);
            }
            set_color(rend, 220, 30, 30, 255);
            draw_text(rend, font_sm, "REC", SIM_W - 44, 6, 220, 30, 30);
            char rfbuf[32];
            snprintf(rfbuf, sizeof(rfbuf), "%ds",
                     rec_frames / 60);
            draw_text(rend, font_sm, rfbuf, SIM_W - 44, 20, 180, 30, 30);
        }

        /* info panel */
        draw_panel(rend, font_lg, font_sm, &sim, &cam, cam.selected);

        SDL_RenderPresent(rend);

        /* ── capture frame to ffmpeg if recording ──────────────── */
        if(recording && ffmpeg_pipe && pixel_buf){
            SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_BGRA32,
                                 pixel_buf, WIN_W * 4);
            fwrite(pixel_buf, 4, WIN_W * WIN_H, ffmpeg_pipe);
            rec_frames++;
        }
    }

    if(recording && ffmpeg_pipe){
        pclose(ffmpeg_pipe); ffmpeg_pipe = NULL;
    }
    free(pixel_buf);

    TTF_CloseFont(font_lg);
    TTF_CloseFont(font_sm);
    TTF_Quit();
    SDL_DestroyRenderer(rend);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
