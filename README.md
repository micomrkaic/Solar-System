# 3-D Solar System Simulation

A real-time interactive simulation of the solar system written in C17 using SDL2.
All eight planets and the Sun are integrated as a full N-body gravitational
system with a 4th-order Yoshida symplectic integrator, with initial conditions
from JPL Horizons (J2000.0 epoch). One source file builds a native binary and a
WebAssembly version; a GitHub Actions workflow runs the physics tests and
deploys the web build to GitHub Pages on every push.

---

## Features

- **4th-order symplectic N-body physics** — Yoshida's composition of three
  kick-drift-kick leapfrog steps (dt = 2 h). Energy is conserved to ~10⁻¹²
  over a century; Mercury's *numerical* perihelion drift is ~1.5″/century.
  Plain 2nd-order leapfrog is selectable (**I**) for comparison — its drift is
  ~22,000″/century at the same step.
- **General relativity** (**G**) — first post-Newtonian Sun–planet term. The
  panel shows Mercury's measured secular perihelion rate ϖ̇ ± σ: ≈ 530″/century
  Newtonian, ≈ 573″/century with GR, the difference being Einstein's 43″.
- **Conservation readouts** — ΔE/E₀ and ΔL/L₀ since t = 0, live.
- **Calendar date** (TT) from the J2000 epoch, forwards or backwards.
- **Time reversal** (**T**) and **pause** (**Space**) — the integrator is
  time-symmetric, so the system retraces its path.
- **Test particles** (**A**) — 1500 main-belt and 500 Kuiper-belt massless
  bodies, integrated under the planets (O(N·M)); Kirkwood-style structure
  develops over time.
- **Initial conditions** — exact JPL Horizons state vectors at J2000.0
  (`tests/fetch_j2000_state.sh` generates `horizons_j2000.h`, picked up by
  `make` automatically). Without the header, the JPL SSD *mean* elements
  table is used; those are not osculating elements and give %-level
  position errors for Saturn after 20 years (Jupiter–Saturn great
  inequality), which the Horizons test then reports rather than fails.
- **Full mutual perturbations** — all bodies exert gravity on all others.
- **Crisp at any size** — the scene is rendered at the display's physical
  resolution (resizable, HiDPI, any browser zoom); nothing is resampled.
- **3-D projection** with perspective, arcball rotation, pan, zoom, and a
  **follow mode** (**F**) that centres a body and fits its orbit.
- **Orbit ellipses**, **fading trails** (one orbital period each), **axis
  snap views**, **click-to-select**, **label de-confliction**, **help
  overlay** (**H**).
- **Video recording** (native) — press **V** to pipe frames into ffmpeg.
- **Self-tests** (`make test`) — vis-viva, Kepler III, conservation, Mercury's
  1PN advance, and optionally positions against JPL Horizons.

---

## Dependencies

| Library | Purpose |
|---------|---------|
| SDL2 | window, renderer, events |
| SDL2_ttf | text rendering |
| ffmpeg | video recording (runtime only, not a compile dependency) |

### Installing dependencies

**macOS (Homebrew)**
```bash
brew install sdl2 sdl2_ttf ffmpeg
```

**Debian / Ubuntu**
```bash
sudo apt install libsdl2-dev libsdl2-ttf-dev ffmpeg
```

**Arch Linux**
```bash
sudo pacman -S sdl2 sdl2_ttf ffmpeg
```

**Fedora**
```bash
sudo dnf install SDL2-devel SDL2_ttf-devel ffmpeg
```

---

## Building

One source file builds both the native binary and the browser version.

```bash
make            # native binary ./solar_system
make run        # build and launch
make info       # show detected build settings
```

The Makefile auto-detects macOS (Homebrew) vs Linux. Manual compile:

```bash
gcc -O2 -std=c17 solar_system.c -o solar_system \
    $(pkg-config --cflags --libs sdl2 SDL2_ttf) -lm
```

### Installing

```bash
make install            # → ~/.local/bin, ~/.local/share/solar_system, ~/.local/share/man/man6
make install PREFIX=/usr/local   # system-wide (needs sudo)
make uninstall
```

This installs the binary, the bundled DejaVu font (so the installed program
needs no system fonts) and a man page: `man solar_system`.

### Web (WebAssembly)

Requires the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html):

```bash
source ~/emsdk/emsdk_env.sh
make wasm       # → web/index.html  (one self-contained file: JS + Wasm + font)
make serve      # test at http://localhost:8080/
```

**GitHub Pages** — two ways:

1. *Recommended:* set *Settings → Pages → Source* to **GitHub Actions**. The
   included workflow builds and deploys on every push to `main`, after the
   physics tests pass. Nothing to commit.
2. *From a branch:* `make pages` copies `web/index.html` → `docs/index.html`;
   commit it and set *Settings → Pages → Source: Deploy from a branch → main,
   /docs*.

The page is at `https://micomrkaic.github.io/solar-system/`.

**itch.io** — `make zip` produces `solar_system_itchio.zip`; upload it as an
HTML project with *This file will be played in the browser* ticked and a
viewport of 1280×720.

The build is a single HTML file on purpose: GitHub Pages serves `.wasm` and
`.data` side files with the wrong MIME type, and a single file is a one-item
upload for itch.io.

---

## Running

```bash
./solar_system
# or
make run
```

The simulation starts at J2000.0 (1 January 2000, 12:00 TT) with all planets
at their true positions and velocities. The native window is resizable; the
scene scales to fit.

The window is resizable and HiDPI-aware; the scene is re-rendered at the new
physical resolution (letterboxed if the aspect ratio differs).

Environment variables: `SOLAR_FONT=/path/to/font.ttf` overrides the font;
`SOLAR_SCALE=1.5` opens a larger window; `SOLAR_SNAPSHOT=frame.bmp`
(optionally with `SOLAR_SNAPSHOT_KEYS=az` and `SOLAR_SNAPSHOT_FRAMES=100`)
renders headlessly, saves the frame and exits — combine with
`SDL_VIDEODRIVER=dummy` on a machine without a display.

---

## Controls

| Input | Action |
|-------|--------|
| **Drag** | Rotate view (arcball) |
| **Shift + drag** | Pan (turns follow mode off) |
| **Scroll wheel** | Zoom in / out |
| **Click a body** | Select it (clicking a name in the panel also works) |
| **1** – **9** | Select body (Sun through Neptune) |
| **F** | Follow mode: centre the selected body, fit its orbit, track it |
| **Z** / **X** / **Y** | Top-down (XY) / edge-on (YZ) / edge-on (XZ) view |
| **R** | Reset view, zoom and pan |
| **Space** | Pause / resume |
| **T** | Reverse time |
| **+** / **-** | Faster / slower (integration steps per frame) |
| **G** | General relativity (1PN Sun–planet term) on / off |
| **I** | Integrator: Yoshida-4 ⇄ Leapfrog-2 |
| **A** | Asteroid + Kuiper belt test particles on / off |
| **H** | Help overlay |
| **V** | Start / stop video recording (native) |
| **Q** or **Esc** | Quit |

---

## Simulation speed

The time step is 2 hours. At startup the simulation runs 10 steps per rendered
frame = 0.83 simulated days per frame (~50 days/second at 60 fps). **+** and
**-** change the steps per frame (1 … 6000, i.e. up to 500 days per frame).
The panel header shows the resulting speed in simulated days or years per real
second, computed from the measured frame rate.

---

## Video recording

Native build only (the browser cannot pipe to ffmpeg). Press **V** to begin recording. A blinking **● REC** indicator appears in the
top-right corner of the simulation area. Press **V** again to stop. The output
file is written to the current directory as:

```
solar_YYYYMMDD_HHMMSS.mp4
```

Encoding: H.264, CRF 18 (near-lossless), 60 fps, yuv420p. Requires `ffmpeg`
in your PATH. `SDL_RenderReadPixels` is used to grab each frame, which adds a
small GPU readback cost; a slight framerate drop while recording is normal.

---

## Physics notes

### Integrator

The basic step is the Störmer-Verlet leapfrog in kick-drift-kick form (velocity
Verlet):

```
v(t + dt/2) = v(t)      + a(t)      · dt/2   [kick]
x(t + dt)   = x(t)      + v(t+dt/2) · dt     [drift]
v(t + dt)   = v(t+dt/2) + a(t+dt)   · dt/2   [kick]
```

The default integrator is Yoshida's 4th-order composition of three such steps
with lengths w₁·dt, w₀·dt, w₁·dt, where w₁ = 1/(2−2^{1/3}) and w₀ = 1−2w₁.
Both are symplectic: they exactly preserve a modified Hamiltonian close to the
true one, so energy oscillates around a conserved value rather than drifting.
The 4th-order scheme matters for the Mercury experiment: at dt = 2 h its
numerical perihelion drift is ~1.5″/century, while plain leapfrog's is
~22,000″/century — 500× the general-relativistic signal it would otherwise
mask. Press **I** to watch the readout switch between the two.

Velocities are kept synchronised with positions, as velocity Verlet requires.
(An earlier version applied the half-step-backwards velocity kick that belongs
to the *staggered* leapfrog; with KDK that biases every planet's eccentricity
vector by ~a·dt/2 — for Mercury, 2×10⁻⁴ in e and 2.7° of perihelion within ten
orbits. `make test` would now catch this.)

### Relativity

With **G** on, each planet receives the dominant first post-Newtonian
acceleration in the Sun's field,

```
a₁PN = μ/(c²r³) · [ (4μ/r − v²) r + 4 (r·v) v ],   μ = G(M☉ + m)
```

which produces Mercury's extra 42.98″/century. The panel measures the secular
rate as a least-squares slope of the osculating longitude of perihelion
(from the Laplace–Runge–Lenz vector) against time; it needs a few decades of
simulated time to settle to ±1″/century. Planet–planet 1PN terms are omitted.

### Units

| Quantity | Unit |
|----------|------|
| Length | AU (astronomical units) |
| Mass | M☉ (solar masses) |
| Time | years |
| G | 4π² AU³ M☉⁻¹ yr⁻² |
| c | 63241.077 AU/yr |

With G = 4π², a circular orbit at 1 AU has period exactly 1 year by
Kepler's third law, which is a useful sanity check.

### Initial conditions

Orbital elements are from the JPL SSD "Keplerian Elements for Approximate
Positions of the Major Planets" (Table 1, J2000.0 epoch):

- Semi-major axis a, eccentricity e, inclination i — taken directly
- Argument of perihelion ω = ϖ − Ω (longitude of perihelion minus longitude
  of ascending node) — this distinction matters; mixing ϖ and ω corrupts
  the orbit orientation by Ω
- Mean anomaly M₀ = L − ϖ (mean longitude minus longitude of perihelion)

Initial velocities in the orbital plane are computed as:

```
vfac = n · a² / r        (n = √(μ/a³), μ = G(M☉+m), r = current radius)
vx   = −vfac · sin(E)
vy   =  vfac · √(1−e²) · cos(E)
```

The factor `n·a²/r` is the correct Keplerian expression. A common error is to
use `n·a/√(1−e²)` instead, which equals `n·a²/r` only at the semi-latus rectum
(r = a(1−e²)) and is wrong everywhere else, giving up to ~23% velocity error for
Mercury at aphelion.

Using n = √(μ/a³) rather than the catalogue period makes the initial state
consistent with the dynamics actually integrated: the vis-viva test passes to
5×10⁻¹⁶.

### Test particles

Asteroids are massless bodies integrated by leapfrog with a 1-day step under
the nine massive bodies; they feel the planets but do not perturb them.
Elements are drawn deterministically (main belt 2.1–3.3 AU, e ≤ 0.15,
i ≤ 12°; Kuiper belt 39–48 AU).

---

## Tests

```bash
make test        # = ./solar_system --test, no window needed
```

| Test | Asserts |
|------|---------|
| vis-viva at t = 0 | \|v² − μ(2/r − 1/a)\|/v² < 10⁻⁹ for every planet |
| Earth two-body period | matches 2π√(a³/μ) to 10⁻⁶ |
| 100-yr N-body conservation | \|ΔE/E₀\| < 10⁻⁸, \|ΔL/L₀\| < 10⁻¹⁰ |
| Mercury numerical drift | \|ϖ̇\| < 5″/century in the two-body problem |
| Mercury 1PN advance | ϖ̇(GR) − ϖ̇(Newton) = 42.98 ± 1″/century |
| Positions vs JPL Horizons | \|Δr\| < 0.1 % of a at the reference epoch (strict with `horizons_j2000.h`; informational from mean elements) |

The Horizons check runs when `tests/horizons.csv` exists;
`tests/fetch_horizons.sh [JD]` downloads heliocentric ecliptic-J2000 positions
of the eight planetary-system barycentres from the Horizons API (default
epoch 2020-01-01). Run `tests/fetch_j2000_state.sh` once and commit the
resulting `horizons_j2000.h` to start from exact osculating states; the
check is then strict.

## Continuous integration

`.github/workflows/ci.yml` builds the native binary, runs `make test`, renders
a headless frame, builds the WebAssembly page and deploys it to GitHub Pages
on every push to `main`. To use it, set *Settings → Pages → Source* to
**GitHub Actions**; `docs/` and `make pages` are then unnecessary (they remain
for deploying from a branch).

---

## File structure

```
solar_system.c           — the whole simulation; native and WebAssembly from one file
Makefile                 — build, test, install/uninstall, wasm, pages, zip
shell.html               — HTML shell the WebAssembly build is wrapped in
man/solar_system.6       — man page (installed by make install)
tests/fetch_horizons.sh  — fetches JPL Horizons reference positions for make test
tests/fetch_j2000_state.sh — generates horizons_j2000.h (exact J2000 initial state)
.github/workflows/ci.yml — tests + Pages deployment
DejaVuSans.ttf           — bundled font (embedded in the web build, installed natively)
docs/index.html          — web build for GitHub Pages when deploying from a branch
```

Platform differences are confined to `#ifdef __EMSCRIPTEN__` blocks in the
source: logical layout size (1600×900 native, 1280×720 web), font location,
the main-loop driver, the resize source, and the ffmpeg recorder.

---

## References

- Bate, Mueller & White, *Fundamentals of Astrodynamics* (1971) — orbital mechanics
- Leimkuhler & Reich, *Simulating Hamiltonian Dynamics* (2004) — symplectic integrators
- JPL Solar System Dynamics, https://ssd.jpl.nasa.gov/planets/approx_pos.html — orbital elements
- Yoshida, "Construction of higher order symplectic integrators", Phys. Lett. A 150 (1990) 262
- Will, *Theory and Experiment in Gravitational Physics* — the 1PN equations of motion
