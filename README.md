# 3-D Solar System Simulation

A real-time interactive simulation of the solar system written in C17 using SDL2.
All eight planets and the Sun are integrated as a full N-body gravitational system
using a symplectic leapfrog (Störmer-Verlet) integrator with initial conditions
from JPL Horizons (J2000.0 epoch).

---

## Features

- **Symplectic N-body physics** — Störmer-Verlet leapfrog integrator preserves
  the symplectic structure of Hamilton's equations, giving excellent long-term
  energy conservation without secular drift.
- **Correct initial conditions** — orbital elements (a, e, i, Ω, ω, M₀) taken
  directly from JPL SSD Table 1 (J2000.0). Velocities computed from the exact
  Keplerian formula v = na²/r, verified against the vis-viva equation.
- **Full mutual perturbations** — all bodies exert gravity on all others; Jupiter's
  effect on Saturn is included.
- **3-D projection** with perspective and arcball (trackball) rotation.
- **Orbit ellipses** drawn analytically from the Keplerian elements, independent
  of the particle positions.
- **Fading trails** per planet, each capped at exactly one orbital period of history.
- **Axis-aligned snap views** — top-down ecliptic (XY), and two edge-on views (XZ, YZ).
- **Live planet data panel** — orbital elements, current distance from Sun, live
  orbital speed, mean orbital speed, radius, and mass.
- **Video recording** — press V to pipe raw frames into ffmpeg and write an MP4.

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

```bash
make
```

The Makefile auto-detects macOS (Homebrew) vs Linux and sets include/library
paths accordingly. A debug build is available with:

```bash
make DEBUG=1
```

To see the detected build settings:

```bash
make info
```

Manual compile if you prefer:

```bash
gcc -O2 -std=c17 solar_system.c -o solar_system \
    $(pkg-config --cflags --libs sdl2 SDL2_ttf) -lm
```

---

## Running

```bash
./solar_system
# or
make run
```

The simulation starts at J2000.0 (1 January 2000, 12:00 TT) with all planets
at their true positions and velocities.

---

## Controls

| Input | Action |
|-------|--------|
| **Drag** | Rotate view (arcball) |
| **Shift + drag** | Pan (translate) |
| **Scroll wheel** | Zoom in / out |
| **+** / **-** | Increase / decrease simulation speed |
| **Z** | Top-down view (XY ecliptic plane) |
| **X** | Edge-on view (YZ plane) |
| **Y** | Edge-on view (XZ plane) |
| **R** | Reset view, zoom, and pan |
| **1** – **9** | Select body (Sun through Neptune) |
| **V** | Start / stop video recording |
| **Q** or **Esc** | Quit |

Clicking a planet name in the right panel also selects it.

---

## Simulation speed

The default timestep is 1/6 day. At startup the simulation runs at 5 steps per
rendered frame, giving roughly 0.83 simulated days per frame (~50 days/second at
60 fps). Use **+** and **-** to adjust; the FPS counter in the bottom-right shows
the current days/frame.

---

## Video recording

Press **V** to begin recording. A blinking **● REC** indicator appears in the
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

The Störmer-Verlet leapfrog in kick-drift-kick form:

```
v(t + dt/2) = v(t)      + a(t)      · dt/2   [kick]
x(t + dt)   = x(t)      + v(t+dt/2) · dt     [drift]
v(t + dt)   = v(t+dt/2) + a(t+dt)   · dt/2   [kick]
```

This is a symplectic integrator: it exactly preserves a modified Hamiltonian
close to the true one, so energy oscillates around a conserved value rather
than drifting. Long-duration orbital stability is qualitatively better than
with Runge-Kutta methods of the same order.

### Units

| Quantity | Unit |
|----------|------|
| Length | AU (astronomical units) |
| Mass | M☉ (solar masses) |
| Time | years |
| G | 4π² AU³ M☉⁻¹ yr⁻² |

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
vfac = n · a² / r        (n = 2π/T, r = current radius)
vx   = −vfac · sin(E)
vy   =  vfac · √(1−e²) · cos(E)
```

The factor `n·a²/r` is the correct Keplerian expression. A common error is to
use `n·a/√(1−e²)` instead, which equals `n·a²/r` only at the semi-latus rectum
(r = a(1−e²)) and is wrong everywhere else, giving up to ~23% velocity error for
Mercury at aphelion.

The leapfrog is initialized with a half-step velocity kick backward in time.
All planet positions must be set before this kick is computed, since it requires
the full N-body acceleration — computing it inside the initialization loop would
use a partially populated position array and corrupt every planet's initial
velocity.

---

## File structure

```
solar_system.c   — full simulation (single file, ~900 lines)
Makefile         — cross-platform build (Linux + macOS/Homebrew)
README.md        — this file
```

---

## References

- Bate, Mueller & White, *Fundamentals of Astrodynamics* (1971) — orbital mechanics
- Leimkuhler & Reich, *Simulating Hamiltonian Dynamics* (2004) — symplectic integrators
- JPL Solar System Dynamics, https://ssd.jpl.nasa.gov/planets/approx_pos.html — orbital elements
