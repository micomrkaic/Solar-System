#define _GNU_SOURCE
/*
 * solar_system_wasm.c  –  3-D Solar System Simulation  (WebAssembly / itch.io port)
 *
 * Ported from solar_system.c for Emscripten.
 * Changes from the native build:
 *   - Window size reduced to 1280×720 for browser embedding.
 *   - Per-frame state moved into a global `App` struct so it survives
 *     across Emscripten's callback-driven main loop.
 *   - Loop body factored into loop_iter(); controlled via
 *     emscripten_set_main_loop() in __EMSCRIPTEN__ builds.
 *   - Font loaded from Emscripten's virtual FS (preloaded DejaVuSans.ttf).
 *   - Video recording (ffmpeg pipe) compiled out in WASM builds.
 *
 * Build: see Makefile.emscripten
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __EMSCRIPTEN__
#  include <emscripten/emscripten.h>
#endif

/* ── constants ─────────────────────────────────────────────────────────── */
#define G_AU        (4.0 * M_PI * M_PI)   /* AU³/(M☉·yr²)               */
#define N_BODIES    9                       /* Sun + 8 planets             */
#define TRAIL_LEN   21000                   /* trail history points        */
#define WIN_W       1280                    /* reduced from 1600 for web   */
#define WIN_H       720                     /* reduced from 900  for web   */
#define PANEL_W     360                     /* right-side info panel width */
#define SIM_W       (WIN_W - PANEL_W)

/* ── body data ─────────────────────────────────────────────────────────── */
typedef struct {
    const char *name;
    double mass;
    double a;
    double ecc;
    double inc;
    double omega;
    double w;
    double M0;
    double radius_km;
    double period_yr;
    double v_orb_km_s;
    Uint8  r, g, b;
    double draw_r;
} BodyDef;

static const BodyDef BODIES[N_BODIES] = {
/*  name         mass(M☉)      a(AU)        ecc         inc(°)     Ω(°)        ω(°)       M0(°)      R(km)      T(yr)   v(km/s) R    G    B   dr */
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
    Vec3 vel[N_BODIES];
    double mass[N_BODIES];

    Vec3  trail[N_BODIES][TRAIL_LEN];
    int   trail_head[N_BODIES];
    int   trail_count[N_BODIES];
    int   trail_cap[N_BODIES];

    double time_yr;
    double dt;
    double speed;
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

static void elements_to_cartesian(const BodyDef *bd, Vec3 *pos, Vec3 *vel){
    if(bd->a == 0.0){ *pos=(Vec3){0,0,0}; *vel=(Vec3){0,0,0}; return; }

    double a   = bd->a;
    double e   = bd->ecc;
    double inc = bd->inc   * M_PI/180.0;
    double Om  = bd->omega * M_PI/180.0;
    double w   = bd->w    * M_PI/180.0;
    double M0  = bd->M0   * M_PI/180.0;

    double E   = solve_kepler(M0, e);

    double nu  = 2.0*atan2(sqrt(1+e)*sin(E/2.0), sqrt(1-e)*cos(E/2.0));
    double r   = a*(1 - e*cos(E));
    double px  = r*cos(nu), py = r*sin(nu);

    double n   = 2.0*M_PI / bd->period_yr;
    double vfac = n*a*a / r;
    double vx  = -vfac*sin(E);
    double vy  =  vfac*sqrt(1-e*e)*cos(E);

    double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);

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
            double d2 = v3dot(r,r) + 1e-10;
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

    for(int i=0;i<N_BODIES;i++)
        sim->vel[i] = v3add(sim->vel[i], v3scale(acc[i], hdt));
    for(int i=0;i<N_BODIES;i++)
        sim->pos[i] = v3add(sim->pos[i], v3scale(sim->vel[i], sim->dt));

    compute_accel(sim->pos, sim->mass, acc);

    for(int i=0;i<N_BODIES;i++)
        sim->vel[i] = v3add(sim->vel[i], v3scale(acc[i], hdt));

    sim->time_yr += sim->dt;
}

/* ── simulation init ─────────────────────────────────────────────────────── */
static void sim_init(Sim *sim){
    memset(sim,0,sizeof(*sim));
    sim->dt    = 1.0/(365.25*6.0);
    sim->speed = 1.0;

    for(int i=0;i<N_BODIES;i++){
        sim->mass[i] = BODIES[i].mass;
        elements_to_cartesian(&BODIES[i], &sim->pos[i], &sim->vel[i]);
    }
    {
        Vec3 acc[N_BODIES];
        compute_accel(sim->pos, sim->mass, acc);
        for(int i=0;i<N_BODIES;i++)
            sim->vel[i] = v3sub(sim->vel[i], v3scale(acc[i], 0.5*sim->dt));
    }
    for(int i=0;i<N_BODIES;i++){
        if(BODIES[i].period_yr <= 0.0){
            sim->trail_cap[i] = 0;
        } else {
            double samples_per_yr = 365.25 / 3.0;
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

/* ── camera / view state ─────────────────────────────────────────────────── */
typedef struct {
    Mat3   rot;
    double scale;
    double fov_dist;
    int    dragging;
    int    panning;
    int    last_x, last_y;
    double pan_x, pan_y;
    int    selected;
} Camera;

static Camera cam_init(void){
    Camera c;
    c.rot     = mat3_identity();
    c.scale   = 55.0;     /* slightly tighter fit for smaller window */
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
    Vec3 ax;
    double angle;
    switch(axis){
        case 0: ax=(Vec3){0,1,0}; angle=-M_PI/2.0; m=mat3_rot_axis(ax,angle); break;
        case 1: ax=(Vec3){1,0,0}; angle= M_PI/2.0; m=mat3_rot_axis(ax,angle); break;
        case 2: m=mat3_identity(); break;
    }
    c->rot = m;
}

static void project(const Camera *cam, Vec3 world,
                    float *sx, float *sy, float *sz_out){
    Vec3 v = mat3_mul_vec(&cam->rot, world);
    double z_ofs = cam->fov_dist;
    double w     = z_ofs / (z_ofs + v.z + 40.0);
    *sx = (float)(SIM_W/2.0 + cam->pan_x + v.x * cam->scale * w);
    *sy = (float)(WIN_H/2.0 + cam->pan_y - v.y * cam->scale * w);
    if(sz_out) *sz_out = (float)v.z;
}

/* ── SDL helpers ─────────────────────────────────────────────────────────── */
static void set_color(SDL_Renderer *rend, Uint8 r, Uint8 g, Uint8 b, Uint8 a){
    SDL_SetRenderDrawColor(rend, r, g, b, a);
}

static void fill_circle(SDL_Renderer *rend, int cx, int cy, int rad){
    for(int dy=-rad; dy<=rad; dy++){
        int dx = (int)sqrt((double)(rad*rad - dy*dy));
        SDL_RenderDrawLine(rend, cx-dx, cy+dy, cx+dx, cy+dy);
    }
}

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
    (void)r; (void)g; (void)b;
    if(bd->a == 0.0) return;
    int STEPS = 180;
    float px=0, py=0;
    for(int step=0; step<=STEPS; step++){
        double nu = 2.0*M_PI*step/STEPS;
        double dist = bd->a*(1.0 - bd->ecc*bd->ecc)/(1.0 + bd->ecc*cos(nu));

        double inc = bd->inc   * M_PI/180.0;
        double Om  = bd->omega * M_PI/180.0;
        double w   = bd->w    * M_PI/180.0;
        double cO=cos(Om),sO=sin(Om),ci=cos(inc),si=sin(inc),cw=cos(w),sw=sin(w);
        double Qxx = cO*cw - sO*sw*ci,  Qxy = -cO*sw - sO*cw*ci;
        double Qyx = sO*cw + cO*sw*ci,  Qyy = -sO*sw + cO*cw*ci;
        double Qzx = sw*si,              Qzy =  cw*si;
        double ox = dist*cos(nu), oy = dist*sin(nu);
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
    (void)cam;
    int px0 = SIM_W + 8;

    set_color(rend,15,15,25,255);
    SDL_Rect bg = {SIM_W, 0, PANEL_W, WIN_H};
    SDL_RenderFillRect(rend, &bg);

    set_color(rend,60,60,80,255);
    SDL_RenderDrawLine(rend, SIM_W, 0, SIM_W, WIN_H);

    draw_text(rend, font_lg, "Solar System", px0, 8, 255,220,50);

    char buf[128];
    double yr = sim->time_yr;
    int years = (int)yr;
    int days  = (int)((yr - years)*365.25);
    snprintf(buf, sizeof(buf), "T = %d yr %d d", years, days);
    draw_text(rend, font_sm, buf, px0, 36, 180,180,200);

    snprintf(buf, sizeof(buf), "Speed: %.0f\xc3\x97 real", sim->speed * sim->dt / (1.0/365.25));
    draw_text(rend, font_sm, buf, px0, 52, 180,180,200);

    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, 72, WIN_W-4, 72);

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

    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, y0+4, WIN_W-4, y0+4);

    if(sel >= 0 && sel < N_BODIES){
        const BodyDef *bd = &BODIES[sel];
        int yd = y0 + 14;
        draw_text(rend, font_lg, bd->name, px0, yd, bd->r, bd->g, bd->b);
        yd += 26;

        Vec3 r = v3sub(sim->pos[sel], sim->pos[0]);
        double dist_AU = v3len(r);

        /* value strings */
        char v_sma[32], v_ecc[32], v_inc[32], v_per[32];
        char v_dist[32], v_spd[32], v_mspd[32], v_rad[32], v_mass[32];

        snprintf(v_sma,  sizeof(v_sma),  "%.4f AU",   bd->a);
        snprintf(v_ecc,  sizeof(v_ecc),  "%.5f",      bd->ecc);
        snprintf(v_inc,  sizeof(v_inc),  "%.3f\xc2\xb0", bd->inc);
        if(bd->period_yr > 0)
            snprintf(v_per, sizeof(v_per), "%.3f yr",  bd->period_yr);
        else
            snprintf(v_per, sizeof(v_per), "\xe2\x80\x93");
        snprintf(v_dist, sizeof(v_dist), "%.4f AU",   dist_AU);

        Vec3 vrel  = v3sub(sim->vel[sel], sim->vel[0]);
        double v_km_s = v3len(vrel) * 4.740;
        snprintf(v_spd,  sizeof(v_spd),  "%.2f km/s", v_km_s);
        snprintf(v_mspd, sizeof(v_mspd), "%.2f km/s", bd->v_orb_km_s);
        snprintf(v_rad,  sizeof(v_rad),  "%.0f km",   bd->radius_km);
        snprintf(v_mass, sizeof(v_mass), "%.3e M\xe2\x8a\x99", bd->mass);

        /* each row: label on left, value at fixed x offset on the same line */
        struct { const char *lbl; const char *val; int section; } rows[] = {
            {"Orbital Elements",   NULL,    1},
            {"  Semi-major axis",  v_sma,   0},
            {"  Eccentricity",     v_ecc,   0},
            {"  Inclination",      v_inc,   0},
            {"  Period",           v_per,   0},
            {"Current State",      NULL,    1},
            {"  Dist. from Sun",   v_dist,  0},
            {"  Speed (now)",      v_spd,   0},
            {"  Speed (mean)",     v_mspd,  0},
            {"Physical",           NULL,    1},
            {"  Radius",           v_rad,   0},
            {"  Mass",             v_mass,  0},
        };

        /* x offset at which values start — leaves room for the longest label */
        int val_x = 128;

        for(int k=0; k<(int)(sizeof(rows)/sizeof(rows[0])); k++){
            if(rows[k].section){
                draw_text(rend, font_sm, rows[k].lbl, px0, yd, 180,180,100);
                yd += 18;
            } else {
                draw_text(rend, font_sm, rows[k].lbl, px0, yd, 140,140,160);
                if(rows[k].val && rows[k].val[0])
                    draw_text(rend, font_sm, rows[k].val, px0+val_x, yd, 220,220,240);
                yd += 16;
            }
        }
    }

    /* controls hint */
    int yh = WIN_H - 100;
    set_color(rend,50,50,70,255);
    SDL_RenderDrawLine(rend, SIM_W+4, yh-4, WIN_W-4, yh-4);
    draw_text(rend, font_sm, "Controls",             px0, yh,    160,160,100);
    draw_text(rend, font_sm, "Drag        \xe2\x80\x93 rotate",  px0, yh+14, 120,120,140);
    draw_text(rend, font_sm, "Shift+drag  \xe2\x80\x93 pan",     px0, yh+28, 120,120,140);
    draw_text(rend, font_sm, "Scroll      \xe2\x80\x93 zoom",    px0, yh+42, 120,120,140);
    draw_text(rend, font_sm, "+/-  \xe2\x80\x93 faster/slower",  px0, yh+56, 120,120,140);
    draw_text(rend, font_sm, "1-9  \xe2\x80\x93 select body",    px0, yh+70, 120,120,140);
    draw_text(rend, font_sm, "R    \xe2\x80\x93 reset view",     px0, yh+84, 120,120,140);
}

/* ── draw legend ───────────────────────────────────────────────────────── */
static void draw_legend(SDL_Renderer *rend, TTF_Font *font,
                        double sim_yr, double speed_mult){
    (void)sim_yr; (void)speed_mult;
    draw_text(rend, font, "3D Solar System \xe2\x80\x93 St\xc3\xb6rmer-Verlet", 8, 6, 100,100,120);
}

/* ════════════════════════════════════════════════════════════════════════
   Application state — lives in globals so loop_iter() can access it
   without needing static locals or a context pointer.
   ════════════════════════════════════════════════════════════════════════ */
static struct {
    SDL_Window   *win;
    SDL_Renderer *rend;
    TTF_Font     *font_lg;
    TTF_Font     *font_sm;
    Camera        cam;
    int           trail_tick;
    int           steps_per_frame;
    int           running;
    Uint32        last_frame;
#ifndef __EMSCRIPTEN__
    FILE  *ffmpeg_pipe;
    Uint8 *pixel_buf;
    int    recording;
    int    rec_frames;
#endif
} app;

/* Sim is large (~6 MB for the trail buffers) – keep it a named static */
static Sim g_sim;

/* ── per-frame callback ─────────────────────────────────────────────────── */
static void loop_iter(void)
{
    SDL_Renderer *rend    = app.rend;
    TTF_Font     *font_lg = app.font_lg;
    TTF_Font     *font_sm = app.font_sm;
    Sim          *sim     = &g_sim;
    Camera       *cam     = &app.cam;

    /* ── events ──────────────────────────────────────────────────────── */
    SDL_Event ev;
    while(SDL_PollEvent(&ev)){
        if(ev.type == SDL_QUIT){
#ifdef __EMSCRIPTEN__
            emscripten_cancel_main_loop();
#else
            app.running = 0;
#endif
        }
        if(ev.type == SDL_KEYDOWN){
            switch(ev.key.keysym.sym){
            case SDLK_ESCAPE:
            case SDLK_q:
#ifdef __EMSCRIPTEN__
                emscripten_cancel_main_loop();
#else
                app.running = 0;
#endif
                break;
            case SDLK_EQUALS: case SDLK_PLUS:
                app.steps_per_frame++;
                if(app.steps_per_frame > 3000) app.steps_per_frame = 3000;
                break;
            case SDLK_MINUS:
                app.steps_per_frame--;
                if(app.steps_per_frame < 1) app.steps_per_frame = 1;
                break;
            case SDLK_r: app.cam = cam_init(); cam = &app.cam; break;
            case SDLK_x: cam_set_view(cam, 0); break;
            case SDLK_y: cam_set_view(cam, 1); break;
            case SDLK_z: cam_set_view(cam, 2); break;
            case SDLK_1: cam->selected=0; break;
            case SDLK_2: cam->selected=1; break;
            case SDLK_3: cam->selected=2; break;
            case SDLK_4: cam->selected=3; break;
            case SDLK_5: cam->selected=4; break;
            case SDLK_6: cam->selected=5; break;
            case SDLK_7: cam->selected=6; break;
            case SDLK_8: cam->selected=7; break;
            case SDLK_9: cam->selected=8; break;
#ifndef __EMSCRIPTEN__
            case SDLK_v: {
                if(!app.recording){
                    time_t t = time(NULL);
                    struct tm *tm_ = localtime(&t);
                    char fname[64];
                    strftime(fname, sizeof(fname), "solar_%Y%m%d_%H%M%S.mp4", tm_);
                    char cmd[256];
                    snprintf(cmd, sizeof(cmd),
                        "ffmpeg -y -f rawvideo -pixel_format bgra "
                        "-video_size %dx%d -framerate 60 -i pipe:0 "
                        "-c:v libx264 -preset fast -crf 18 "
                        "-pix_fmt yuv420p %s",
                        WIN_W, WIN_H, fname);
                    app.ffmpeg_pipe = popen(cmd, "w");
                    if(app.ffmpeg_pipe){
                        app.pixel_buf = malloc(WIN_W * WIN_H * 4);
                        app.recording = 1; app.rec_frames = 0;
                        SDL_Log("Recording started: %s", fname);
                    }
                } else {
                    app.recording = 0;
                    free(app.pixel_buf); app.pixel_buf = NULL;
                    pclose(app.ffmpeg_pipe); app.ffmpeg_pipe = NULL;
                    SDL_Log("Recording stopped (%d frames)", app.rec_frames);
                }
            } break;
#endif /* !__EMSCRIPTEN__ */
            }
        }
        if(ev.type == SDL_MOUSEBUTTONDOWN && ev.button.button==SDL_BUTTON_LEFT){
            SDL_Keymod mod = SDL_GetModState();
            cam->panning  = (mod & KMOD_SHIFT) ? 1 : 0;
            cam->dragging = 1;
            cam->last_x   = ev.button.x;
            cam->last_y   = ev.button.y;
            /* click in panel selects body */
            int mx=ev.button.x, my=ev.button.y;
            if(mx >= SIM_W){
                int idx = (my - 80) / 22;
                if(idx>=0 && idx<N_BODIES) cam->selected=idx;
            }
        }
        if(ev.type == SDL_MOUSEBUTTONUP && ev.button.button==SDL_BUTTON_LEFT){
            cam->dragging = 0;
            cam->panning  = 0;
        }
        if(ev.type == SDL_MOUSEMOTION && cam->dragging){
            int dx = ev.motion.x - cam->last_x;
            int dy = ev.motion.y - cam->last_y;
            cam->last_x = ev.motion.x; cam->last_y = ev.motion.y;
            if(cam->panning){
                cam->pan_x += dx;
                cam->pan_y += dy;
            } else {
                Vec3 yax={0,1,0}, xax={1,0,0};
                Mat3 Ry = mat3_rot_axis(yax, dx*0.005);
                Mat3 Rx = mat3_rot_axis(xax, dy*0.005);
                cam->rot = mat3_mul(&Ry, &cam->rot);
                cam->rot = mat3_mul(&Rx, &cam->rot);
            }
        }
        if(ev.type == SDL_MOUSEWHEEL){
            cam->scale *= (ev.wheel.y > 0) ? 1.12 : 0.89;
            if(cam->scale < 3.0)    cam->scale = 3.0;
            if(cam->scale > 8000.0) cam->scale = 8000.0;
        }
    }

    /* ── physics ─────────────────────────────────────────────────── */
    for(int step=0; step<app.steps_per_frame; step++){
        leapfrog_step(sim);
        app.trail_tick++;
        if(app.trail_tick >= 18){
            app.trail_tick = 0;
            for(int i=0;i<N_BODIES;i++) trail_push(sim,i);
        }
    }

    /* ── render ──────────────────────────────────────────────────── */
    set_color(rend,5,5,12,255);
    SDL_RenderClear(rend);

    /* starfield */
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

    /* orbit ellipses */
    for(int i=1;i<N_BODIES;i++){
        Uint8 r=BODIES[i].r, g=BODIES[i].g, b=BODIES[i].b;
        SDL_SetRenderDrawColor(rend, r/3, g/3, b/3, 80);
        draw_orbit(rend, cam, &BODIES[i], r/3, g/3, b/3);
    }

    /* trails */
    for(int i=1;i<N_BODIES;i++){
        int cnt = sim->trail_count[i];
        int head = sim->trail_head[i];
        int cap  = sim->trail_cap[i];
        if(cap <= 0) continue;
        float px2=0,py2=0;
        for(int k=0;k<cnt;k++){
            int idx = ((head - cnt + k) + cap) % cap;
            float sx,sy;
            project(cam, sim->trail[i][idx], &sx, &sy, NULL);
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
        project(cam, sim->pos[i], &sx, &sy, NULL);
        int rad = (int)BODIES[i].draw_r;

        if(i==0){
            for(int glow=rad+8; glow>rad; glow--){
                Uint8 a=(Uint8)(30*(glow-rad));
                SDL_SetRenderDrawColor(rend,255,200,50,a);
                for(int dy2=-glow;dy2<=glow;dy2++){
                    int dx2=(int)sqrt((double)(glow*glow-dy2*dy2));
                    SDL_RenderDrawLine(rend,(int)sx-dx2,(int)sy+dy2,(int)sx+dx2,(int)sy+dy2);
                }
            }
        }

        set_color(rend, BODIES[i].r, BODIES[i].g, BODIES[i].b, 255);
        fill_circle(rend, (int)sx, (int)sy, rad);

        if(i == cam->selected){
            SDL_SetRenderDrawColor(rend,255,255,255,160);
            int rr=rad+4;
            for(int ang=0;ang<360;ang+=4){
                float a1=ang*M_PI/180.0f, a2=(ang+4)*M_PI/180.0f;
                SDL_RenderDrawLineF(rend,
                    sx+rr*cosf(a1), sy+rr*sinf(a1),
                    sx+rr*cosf(a2), sy+rr*sinf(a2));
            }
        }

        draw_text(rend, font_sm, BODIES[i].name,
                  (int)sx+rad+4, (int)sy-7,
                  BODIES[i].r, BODIES[i].g, BODIES[i].b);
    }

    /* axis indicator */
    {
        int ax0=50, ay0=WIN_H-50;
        float len=30.0f;
        Vec3 ox={1,0,0},oy={0,1,0},oz={0,0,1};
        Vec3 rx=mat3_mul_vec(&cam->rot,ox);
        Vec3 ry=mat3_mul_vec(&cam->rot,oy);
        Vec3 rz=mat3_mul_vec(&cam->rot,oz);
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

    draw_legend(rend, font_sm, sim->time_yr, (double)app.steps_per_frame);

    /* FPS counter */
    {
        Uint32 now=SDL_GetTicks();
        double fps = 1000.0/(now - app.last_frame + 1);
        app.last_frame = now;
        double days_per_frame = app.steps_per_frame * sim->dt * 365.25;
        char fbuf[64];
        snprintf(fbuf,64,"%.0f fps  %.2f days/frame", fps, days_per_frame);
        draw_text(rend,font_sm,fbuf, SIM_W-200, WIN_H-20, 80,80,100);
    }

#ifndef __EMSCRIPTEN__
    /* REC indicator (native only) */
    if(app.recording){
        Uint32 ms = SDL_GetTicks();
        if((ms/500)%2 == 0){
            set_color(rend, 220, 30, 30, 255);
            fill_circle(rend, SIM_W - 18, 14, 7);
        }
        set_color(rend, 220, 30, 30, 255);
        draw_text(rend, font_sm, "REC", SIM_W - 44, 6, 220, 30, 30);
        char rfbuf[32];
        snprintf(rfbuf, sizeof(rfbuf), "%ds", app.rec_frames / 60);
        draw_text(rend, font_sm, rfbuf, SIM_W - 44, 20, 180, 30, 30);
    }
#endif

    draw_panel(rend, font_lg, font_sm, sim, cam, cam->selected);

    SDL_RenderPresent(rend);

#ifndef __EMSCRIPTEN__
    if(app.recording && app.ffmpeg_pipe && app.pixel_buf){
        SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_BGRA32,
                             app.pixel_buf, WIN_W * 4);
        fwrite(app.pixel_buf, 4, WIN_W * WIN_H, app.ffmpeg_pipe);
        app.rec_frames++;
    }
#endif
}

/* ── main ───────────────────────────────────────────────────────────────── */
int main(void){
    if(SDL_Init(SDL_INIT_VIDEO) < 0){
        fprintf(stderr,"SDL_Init: %s\n", SDL_GetError()); return 1;
    }
    if(TTF_Init() < 0){
        fprintf(stderr,"TTF_Init: %s\n", TTF_GetError()); return 1;
    }

    /* In Emscripten, SDL_WINDOW_RESIZABLE is not set – the canvas size is
       fixed by the HTML shell.  On native builds we keep resizable.       */
    Uint32 win_flags = SDL_WINDOW_SHOWN;
#ifndef __EMSCRIPTEN__
    win_flags |= SDL_WINDOW_RESIZABLE;
#endif

    app.win  = SDL_CreateWindow("Solar System",
                    SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                    WIN_W, WIN_H, win_flags);
    app.rend = SDL_CreateRenderer(app.win, -1,
                    SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if(!app.win || !app.rend){
        fprintf(stderr,"Window/Renderer: %s\n", SDL_GetError()); return 1;
    }

    /* ── font loading ─────────────────────────────────────────────────
       In WASM builds only the preloaded virtual-FS path is tried.
       In native builds a list of system font paths is searched.         */
#ifdef __EMSCRIPTEN__
    const char *font_paths[] = {
        "DejaVuSans.ttf",  /* preloaded via --preload-file               */
        NULL
    };
#else
    const char *font_paths[] = {
        /* macOS */
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
        "/opt/homebrew/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/dejavu/DejaVuSans.ttf",
        /* Linux */
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
        NULL
    };
#endif

    for(int i=0; font_paths[i]; i++){
        app.font_lg = TTF_OpenFont(font_paths[i], 15);
        app.font_sm = TTF_OpenFont(font_paths[i], 12);
        if(app.font_lg && app.font_sm) break;
        if(app.font_lg){ TTF_CloseFont(app.font_lg); app.font_lg=NULL; }
        if(app.font_sm){ TTF_CloseFont(app.font_sm); app.font_sm=NULL; }
    }
    if(!app.font_lg || !app.font_sm){
        fprintf(stderr,
            "Could not open font.\n"
#ifdef __EMSCRIPTEN__
            "Ensure DejaVuSans.ttf was preloaded (--preload-file DejaVuSans.ttf).\n"
#else
            "Install DejaVu fonts or set SOLAR_FONT=/path/to/font.ttf\n"
#endif
        );
#ifndef __EMSCRIPTEN__
        const char *env_font = getenv("SOLAR_FONT");
        if(env_font){
            app.font_lg = TTF_OpenFont(env_font, 15);
            app.font_sm = TTF_OpenFont(env_font, 12);
        }
        if(!app.font_lg || !app.font_sm)
#endif
        return 1;
    }

    /* ── init simulation & camera ────────────────────────────────────── */
    sim_init(&g_sim);
    app.cam             = cam_init();
    app.trail_tick      = 0;
    app.steps_per_frame = 5;
    app.running         = 1;
    app.last_frame      = SDL_GetTicks();

#ifndef __EMSCRIPTEN__
    app.ffmpeg_pipe = NULL;
    app.pixel_buf   = NULL;
    app.recording   = 0;
    app.rec_frames  = 0;
#endif

    SDL_SetRenderDrawBlendMode(app.rend, SDL_BLENDMODE_BLEND);

    /* ── start main loop ─────────────────────────────────────────────── */
#ifdef __EMSCRIPTEN__
    /* 0 fps = browser controls the frame rate (requestAnimationFrame).
       1 = simulate_infinite_loop: let Emscripten manage the loop.       */
    emscripten_set_main_loop(loop_iter, 0, 1);
    /* Code below this point is unreachable in WASM builds. */
#else
    while(app.running){
        loop_iter();
    }

    /* ── cleanup (native only) ─────────────────────────────────────── */
    if(app.recording && app.ffmpeg_pipe){
        pclose(app.ffmpeg_pipe); app.ffmpeg_pipe = NULL;
    }
    free(app.pixel_buf);
    TTF_CloseFont(app.font_lg);
    TTF_CloseFont(app.font_sm);
    TTF_Quit();
    SDL_DestroyRenderer(app.rend);
    SDL_DestroyWindow(app.win);
    SDL_Quit();
#endif
    return 0;
}
