# Makefile — 3-D Solar System Simulation
# Linux (apt/pacman/dnf) and macOS (Homebrew, arm64 + x86_64); Emscripten for the web.
#
# Targets:
#   make             → native binary   ./solar_system
#   make run         → build + launch  (native)
#   make test        → physics self-tests (vis-viva, Kepler III, conservation,
#                      Mercury 1PN precession, optional JPL Horizons check)
#   tests/fetch_j2000_state.sh → horizons_j2000.h: exact JPL initial conditions
#                      (picked up automatically by make; commit the header)
#   make install     → install binary, font and man page to ~/.local (or PREFIX=…)
#   make uninstall   → remove installed files
#   make wasm        → WebAssembly build  web/index.html  (single self-contained file)
#   make serve       → local test server for the WASM build
#   make pages       → copy WASM build to docs/ (GitHub Pages: main branch, /docs)
#   make zip         → solar_system_itchio.zip for itch.io
#   make clean       → remove binaries, web/, zip
#   make info        → show detected build settings
#
# One source file (solar_system.c) builds both targets; see its header.

CC      ?= gcc
EMCC    ?= emcc
PYTHON  ?= python3

TARGET   = solar_system
SRC      = solar_system.c
FONT     = DejaVuSans.ttf
WEB_DIR  = web

# ── install locations ─────────────────────────────────────────────────────
PREFIX  ?= $(HOME)/.local
BINDIR   = $(PREFIX)/bin
DATADIR  = $(PREFIX)/share/solar_system
MANDIR   = $(PREFIX)/share/man/man6

# ── platform detection ────────────────────────────────────────────────────
UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
    BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /usr/local)
    PKG_CONFIG  := PKG_CONFIG_PATH=$(BREW_PREFIX)/lib/pkgconfig pkg-config
else
    PKG_CONFIG  := pkg-config
endif

# ── initial conditions ────────────────────────────────────────────────────
# If tests/fetch_j2000_state.sh has produced horizons_j2000.h, start from
# JPL Horizons' exact J2000 state vectors instead of the mean-elements table.
IC_FLAGS := $(if $(wildcard horizons_j2000.h),-DHAVE_HORIZONS_J2000,)

# ── native flags ──────────────────────────────────────────────────────────
# SOLAR_DATADIR is baked in so an installed binary finds the installed font.
CFLAGS  := -O2 -std=c17 -Wall -Wextra -Wpedantic -Wno-unused-parameter \
           -DSOLAR_DATADIR='"$(DATADIR)"' $(IC_FLAGS) \
           $(shell $(PKG_CONFIG) --cflags sdl2 SDL2_ttf)
LDFLAGS := $(shell $(PKG_CONFIG) --libs sdl2 SDL2_ttf) -lm
ifeq ($(UNAME), Darwin)
    LDFLAGS += -framework Cocoa
endif

# ── WebAssembly flags ─────────────────────────────────────────────────────
# -std=gnu17            EM_ASM-style macros need GNU extensions
# SINGLE_FILE=1         wasm + embedded font baked into index.html as base64:
#                       no .wasm/.data side files, so no MIME-type problems
#                       on GitHub Pages and a one-file upload for itch.io
# --embed-file          font goes into the binary's virtual FS (no XHR)
# GROWABLE_ARRAYBUFFERS=0  emcc 6.0.2 defaulted this to 1, which backs the
#                       heap with a resizable ArrayBuffer that Chrome's
#                       TextDecoder/WebGL reject (reverted upstream in 6.0.3)
# MINIFY_HTML=0         distro emscripten lacks html-minifier-terser
CFLAGS_WASM := -O2 -std=gnu17 $(IC_FLAGS) \
               -s USE_SDL=2 \
               -s USE_SDL_TTF=2 \
               -s ALLOW_MEMORY_GROWTH=1 \
               -s GROWABLE_ARRAYBUFFERS=0 \
               -s INITIAL_MEMORY=67108864 \
               -s ENVIRONMENT=web \
               -s SINGLE_FILE=1 \
               -s ASSERTIONS=0 \
               -s MINIFY_HTML=0 \
               --embed-file $(FONT) \
               --shell-file shell.html

# ── native ────────────────────────────────────────────────────────────────
.PHONY: all run test install uninstall wasm serve pages zip clean info check-emcc check-font

all: $(TARGET)

$(TARGET): $(SRC) $(wildcard horizons_j2000.h)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

test: $(TARGET)
	./$(TARGET) --test

# ── install / uninstall ───────────────────────────────────────────────────
install: $(TARGET)
	install -Dm755 $(TARGET)             $(BINDIR)/$(TARGET)
	install -Dm644 $(FONT)               $(DATADIR)/$(FONT)
	install -Dm644 man/solar_system.6    $(MANDIR)/solar_system.6
	@echo "Installed binary  : $(BINDIR)/$(TARGET)"
	@echo "Installed font    : $(DATADIR)/$(FONT)"
	@echo "Installed man page: $(MANDIR)/solar_system.6   (man solar_system)"
	@echo "Ensure $(BINDIR) is in PATH and $(PREFIX)/share/man is in MANPATH."

uninstall:
	rm -f  $(BINDIR)/$(TARGET)
	rm -rf $(DATADIR)
	rm -f  $(MANDIR)/solar_system.6
	@echo "Uninstalled."

# ── WebAssembly ───────────────────────────────────────────────────────────
wasm: check-emcc check-font
	mkdir -p $(WEB_DIR)
	$(EMCC) $(SRC) $(CFLAGS_WASM) -o $(WEB_DIR)/index.html
	@echo "Built $(WEB_DIR)/index.html — one self-contained file."

serve: wasm
	@echo "Serving at http://localhost:8080/"
	cd $(WEB_DIR) && $(PYTHON) -m http.server 8080

pages: wasm
	mkdir -p docs
	cp $(WEB_DIR)/index.html docs/index.html
	@echo "docs/ is ready — commit it; GitHub Pages: Settings → Pages → main branch, /docs."

zip: wasm
	cd $(WEB_DIR) && zip -j ../solar_system_itchio.zip index.html
	@echo "Created solar_system_itchio.zip — upload to itch.io (kind: HTML)."

# ── housekeeping ──────────────────────────────────────────────────────────
clean:
	$(RM) $(TARGET) solar_system_itchio.zip
	$(RM) -r $(WEB_DIR)

info:
	@echo "Platform : $(UNAME)"
	@echo "Compiler : $(CC)"
	@echo "CFLAGS   : $(CFLAGS)"
	@echo "LDFLAGS  : $(LDFLAGS)"
	@echo "PREFIX   : $(PREFIX)"

# ── guards ────────────────────────────────────────────────────────────────
check-emcc:
	@command -v $(EMCC) >/dev/null 2>&1 || { \
	  echo "ERROR: emcc not found.  Activate the Emscripten SDK first:"; \
	  echo "  source ~/emsdk/emsdk_env.sh"; \
	  exit 1; }

# The font is tracked in the repo; this only matters for a stripped checkout.
check-font:
	@if [ ! -f $(FONT) ]; then \
	  echo "$(FONT) not found locally — searching system font paths…"; \
	  FOUND=""; \
	  for p in \
	    /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf \
	    /usr/share/fonts/dejavu/DejaVuSans.ttf \
	    /usr/share/fonts/TTF/DejaVuSans.ttf \
	    /usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf \
	    /usr/local/share/fonts/dejavu/DejaVuSans.ttf \
	    "$$(brew --prefix 2>/dev/null)/share/fonts/dejavu/DejaVuSans.ttf" \
	    /opt/homebrew/share/fonts/dejavu/DejaVuSans.ttf \
	    /Library/Fonts/DejaVuSans.ttf \
	  ; do \
	    if [ -f "$$p" ]; then FOUND="$$p"; break; fi; \
	  done; \
	  if [ -n "$$FOUND" ]; then echo "Found at $$FOUND — copying…"; cp "$$FOUND" $(FONT); \
	  else \
	    echo "ERROR: Could not find DejaVuSans.ttf."; \
	    echo "  Debian/Ubuntu : sudo apt install fonts-dejavu-core"; \
	    echo "  macOS         : brew install --cask font-dejavu"; \
	    exit 1; \
	  fi; \
	fi
