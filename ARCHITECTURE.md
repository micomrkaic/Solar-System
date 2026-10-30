# Architecture and design decisions

This file records *why* the simulation is built the way it is — the
measurements that drove each decision, so that nobody has to rediscover them.
For how to build, run and use it, see `README.md`; for the controls and physics
summary, `man solar_system`.

## Structure

One C17 source, `solar_system.c`, builds every target. Reading order:

| Section | What it holds |
|---|---|
| constants, `BODIES[]` | units, logical layout size, DE441 masses, JPL mean elements (display + fallback) |
| physics | `kepler_to_cartesian`, `compute_accel` (Newton + 1PN), `kdk_step`, `sim_step` (Yoshida-4), `sim_energy`, `osc_varpi`, `sim_date` |
| `Belt` | massless test particles, own leapfrog at 12×dt |
| camera / device transform | `project()` to logical px; `DX/DY/DL` to device px; `line`, `polyline_dev`, `dot`, `fill_circle`, `draw_text` |
| app state + rendering | `draw_orbit`, `draw_trail`, `draw_bodies` (label placement), `draw_panel`, `draw_help` |
| resolution, input, loop | `apply_resolution`, `poll_resize` (web), `handle_key`, `loop_iter` |
| tests | `run_tests` (`--test`), `dump_state` (`--dump-state`) |
| `main` | window/renderer/font setup, Emscripten vs native loop |

Platform differences live only in `#ifdef __EMSCRIPTEN__` blocks: logical
window size, font location, main-loop driver, resize source, ffmpeg recorder.

## Physics

### Integrator: Yoshida-4 at dt = 2 h

The step is three kick-drift-kick (velocity Verlet) sub-steps of length
w₁·dt, w₀·dt, w₁·dt with w₁ = 1/(2−2^{1/3}), w₀ = 1−2w₁ (Yoshida 1990).

Measured on the Sun–Mercury two-body problem over 100 years, perihelion drift
that should be exactly zero:

| integrator | dt | numerical drift | max \|ΔE/E₀\| |
|---|---|---|---|
| leapfrog-2 | 4 h | −21,928″/century | 3.9e-5 |
| Yoshida-4 | 4 h | +32″/century | 1.4e-8 |
| Yoshida-4 | 2 h | +1.45″/century | 9.6e-11 |
| Yoshida-4 | 1 h | −0.46″/century | 8.9e-10 |

The general-relativistic signal the program is meant to display is
42.98″/century. Plain leapfrog at a usable step buries it 500× over; Yoshida-4
at 2 h is 30× under it at ~6× leapfrog's cost per unit time, which is
negligible for nine bodies. Leapfrog remains selectable (key **I**) precisely
so the −22,000 can be seen.

### Velocity initialisation

KDK keeps velocities synchronised with positions. An earlier version applied
the half-step-backwards kick that belongs to the *staggered* leapfrog; with
KDK that is an O(dt) velocity bias of a·dt/2 — for Mercury, 2×10⁻⁴ in
eccentricity and 2.7° of perihelion longitude within ten orbits. Removed.

### Units and GM☉

AU, solar masses, Julian years (365.25 d). GM☉ = k²·365.25² =
39.476926421373 AU³/yr² with k = 0.01720209895 the Gaussian gravitational
constant — this equals DE441's GM☉ with the IAU 2012 AU.

The textbook G = 4π² is exact only in Gaussian years (365.2569 d) and is
3.8×10⁻⁵ too large here. Because the Horizons initial *velocities* are fixed
numbers, a too-large μ also shrinks the osculating semi-major axis
(δa/a ≈ −δμ/μ), so the mean-motion error is δn/n = ½·δμ/μ + 1.5·δμ/μ =
7.6×10⁻⁵ — ~0.01 AU of along-track error for every inner planet after 20
years. The Horizons position test exposed it as a uniform δn/n across Venus,
Earth, Mars, Jupiter and Saturn.

Masses are DE441 GM ratios of the planetary *systems* (Earth = Earth–Moon
barycentre), consistent with the barycentric Horizons vectors.

### Initial conditions

Horizons state vectors at J2000.0 for bodies 10 and 1–8 (system
barycentres), solar-system-barycentric, ecliptic/J2000, fetched by
`tests/fetch_j2000_state.sh` into `horizons_j2000.h` and committed. The JPL
approximate elements remain as a fallback and for the data panel, but they are
*mean* elements: used as osculating ones they put Saturn 3.5 % off after 20
years through the Jupiter–Saturn great inequality.

### Validation against DE441

With exact initial conditions, GM☉ as above, and the 1PN term on, heliocentric
positions at JD 2458849.5 (2020-01-01) differ from Horizons by:

| Mercury | Venus | Earth | Mars | Jupiter | Saturn | Uranus | Neptune |
|---|---|---|---|---|---|---|---|
| 800 km | 1800 km | 1100 km | 750 km | 113 km | 54 km | 14 km | 15 km |

The remainder is physics the nine-body model omits: ~340 asteroids in DE441,
the Moon as a separate body, solar oblateness. The same omissions are why the
full N-body GR run reads ϖ̇(Mercury) ≈ 569.5″/century instead of the accepted
~574.

### Relativity

First post-Newtonian Sun–planet acceleration
a = μ/(c²r³)·[(4μ/r − v²)r + 4(r·v)v], μ = G(M☉+m), c = 63241.077 AU/yr,
with the reaction on the Sun. Planet–planet 1PN terms (~10⁻⁴ of this) are
omitted. Measured: +42.81″/century on Sun–Mercury (theory 42.98).

The perihelion-rate readout is a least-squares slope of the unwrapped
osculating longitude of perihelion (from the Laplace–Runge–Lenz vector)
against time. Short-period terms from Venus and Jupiter dominate for the first
years, so nothing is shown before 5 simulated years; it settles to
±1″/century after a few decades.

### Test particles

1500 main-belt (2.1–3.3 AU) and 500 Kuiper-belt (39–48 AU) massless bodies,
deterministic elements, integrated by plain leapfrog at 12×dt = 1 day under
the nine massive bodies. Cost O(N·M); at the maximum speed setting with belts
on, frame time rises noticeably — that is expected.

## Rendering

### Device-resolution drawing

The scene is laid out in a fixed logical space (1600×900 native, 1280×720
web). The back-buffer is the display's *physical* pixel size; `g_S` scales
and `g_ox/g_oy` letterbox. Everything is drawn in device pixels:

- fonts are reopened at `size × g_S` so FreeType rasterises at device
  resolution and glyphs are blitted 1:1;
- lines and polylines are drawn in device space with stroke `round(g_S)`,
  offset copies across the minor axis;
- orbits and trails are *polylines* (`SDL_RenderDrawLinesF`), not per-segment
  calls: SDL paints interior joints once, which with alpha blending removes
  the bright dot at every vertex the old code had.

The earlier approach — a fixed 1280×720 canvas downscaled by CSS with
`image-rendering: pixelated` — dropped whole rows and columns and turned 1-px
lines and text strokes into dashes.

Web: `shell.html` measures the largest 16:9 box in device pixels (`_devW`),
C sizes the canvas to it and calls `_fitCss(W,H)` so the CSS size is exactly
W/dpr × H/dpr; resize is polled each frame. Native: SDL_WINDOW_RESIZABLE +
ALLOW_HIGHDPI, `SDL_WINDOWEVENT_SIZE_CHANGED` → `apply_resolution()`.

### Labels

Greedy placement right / below / above of each body, farthest first; a label
that would overlap one already placed is dropped, except for the Sun and the
selected body.

## Toolchain

- **emcc 6.0.2** shipped with `GROWABLE_ARRAYBUFFERS=1` by default, backing the
  Wasm heap with a resizable `ArrayBuffer` that Chrome's `TextDecoder` and
  WebGL reject ("must not be resizable"). Upstream reverted the default in
  6.0.3. We pass `-sGROWABLE_ARRAYBUFFERS=0` *and* make
  `WebAssembly.Memory.prototype.toResizableBuffer` throw in the shell, so the
  runtime's `try { toResizableBuffer() } catch {}` falls back to the fixed
  buffer on any emcc version.
- `SINGLE_FILE=1` + `--embed-file DejaVuSans.ttf`: one HTML file, no `.wasm`
  or `.data` side files, hence no MIME-type trouble on GitHub Pages and a
  one-item upload for itch.io.
- `-std=gnu17` (EM_ASM needs GNU extensions); `MINIFY_HTML=0` (distro
  emscripten lacks html-minifier-terser).

## Testing

`solar_system --test` runs without a window: initial-state sanity, Kepler III
on Sun–Earth, 100-year conservation, Mercury's numerical drift and 1PN
advance, and positions against `tests/horizons.csv` (strict with the J2000
header, informational from mean elements). Headless frames:
`SDL_VIDEODRIVER=dummy SOLAR_SNAPSHOT=f.bmp SOLAR_SNAPSHOT_KEYS=az`.
`--dump-state` prints the current initial conditions in header format.

CI (`.github/workflows/ci.yml`) runs all of this on every push and deploys the
web build to GitHub Pages.
