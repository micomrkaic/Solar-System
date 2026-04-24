# Solar System — WebAssembly / itch.io Build Guide

This document explains how to compile `solar_system_wasm.c` to WebAssembly
and publish the result on itch.io.

---

## Files in this package

| File | Purpose |
|------|---------|
| `solar_system_wasm.c` | Modified C source (Emscripten-ready) |
| `Makefile.emscripten` | Emscripten build rules |
| `index.html` | HTML shell (itch.io iframe host) |
| `solar_system.c` | Original native source (unchanged) |
| `Makefile` | Original native build (unchanged) |

---

## What changed from the native build

| Area | Change |
|------|--------|
| Window size | 1600×900 → **1280×720** (fits itch.io iframe) |
| Main loop | `while(running){}` → **`emscripten_set_main_loop()`** |
| Font loading | System font search → **preloaded `DejaVuSans.ttf`** via Emscripten VFS |
| Video recording | Compiled out with `#ifndef __EMSCRIPTEN__` (no `popen` in WASM) |
| Cleanup code | Skipped in WASM (unreachable after `emscripten_set_main_loop`) |

All physics, rendering, controls, and panel logic are **byte-for-byte identical**
to the native build.

---

## Step 1 — Install the Emscripten SDK

```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh    # re-run in each new terminal session
```

Verify with:
```bash
emcc --version
# expected: emcc (Emscripten ...) 3.x.x
```

---

## Step 2 — Get DejaVuSans.ttf

The font must be in the same directory as `solar_system_wasm.c`.

**Debian / Ubuntu:**
```bash
sudo apt install fonts-dejavu-core
cp /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf .
```

**macOS (Homebrew):**
```bash
brew install --cask font-dejavu
cp "$(brew --prefix)/share/fonts/dejavu/DejaVuSans.ttf" .
```

**Manual download:**
```
https://github.com/dejavu-fonts/dejavu-fonts/releases
```
Download the archive and extract `DejaVuSans.ttf`.

---

## Step 3 — Build

```bash
# Release build (outputs to dist/)
make -f Makefile.emscripten

# Debug build (slower, with ASSERTIONS and SAFE_HEAP)
make -f Makefile.emscripten debug
```

The first run will automatically download SDL2 and SDL2_ttf from Emscripten's
ports system — this takes a few minutes but only happens once.

Build output in `dist/`:
```
dist/
  index.html         ← itch.io HTML shell
  solar_system.js    ← Emscripten JS glue + WASM loader
  solar_system.wasm  ← compiled simulation
  solar_system.data  ← preloaded font file
```

---

## Step 4 — Test locally

Browsers block `file://` origins from loading `.wasm` / `.data`.
Serve with any local HTTP server:

```bash
cd dist
python3 -m http.server 8080
# open http://localhost:8080
```

Or with Node:
```bash
npx serve dist
```

---

## Step 5 — Package and upload to itch.io

```bash
make -f Makefile.emscripten zip
# produces: solar_system_itchio.zip
```

On itch.io:
1. Create a new project → **HTML** kind.
2. Upload `solar_system_itchio.zip`.
3. Tick **"This file will be played in the browser"**.
4. Set viewport to **1280 × 720** (matches `WIN_W` / `WIN_H`).
5. Optionally enable **"Mobile friendly"** (the simulation works with
   touch devices as drag-to-rotate).

---

## Controls (unchanged from native)

| Input | Action |
|-------|--------|
| Drag | Rotate view (arcball) |
| Shift + drag | Pan |
| Scroll wheel | Zoom |
| `+` / `-` | Increase / decrease simulation speed |
| `Z` | Top-down view (XY ecliptic plane) |
| `X` | Edge-on view (YZ plane) |
| `Y` | Edge-on view (XZ plane) |
| `R` | Reset view |
| `1`–`9` | Select body |
| `Q` / `Esc` | Quit (no-op in browser) |

> **Note:** Video recording (`V` key) is disabled in the WASM build because
> it requires `ffmpeg` via `popen()`, which is not available in the browser.

---

## Troubleshooting

**Black screen / nothing loads**
- Check the browser console for errors.
- Make sure you are serving over HTTP, not `file://`.

**"Cannot open font" error in console**
- `DejaVuSans.ttf` was not found in the build directory when compiling.
- Rebuild after copying the font file alongside `solar_system_wasm.c`.

**Very slow / low FPS**
- The N-body integrator runs entirely on the CPU in WASM.
- Reduce simulation speed with `-` or close other browser tabs.
- On mobile/low-power devices, WASM performance is roughly 2–4× slower
  than a native build.

**emcc: command not found**
- Run `source /path/to/emsdk/emsdk_env.sh` to activate the SDK.

---

## Re-building after changes to the C source

Just re-run `make -f Makefile.emscripten`. Emscripten only recompiles
what changed; the SDL ports are already cached.
