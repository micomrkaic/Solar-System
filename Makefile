# Makefile — Solar System Simulator
# Supports: Linux (apt/pacman/dnf) and macOS (Homebrew arm64 + x86_64)
#
# Targets:
#   make          – build
#   make run      – build and launch
#   make clean    – remove binary
#   make info     – show detected build settings

CC     := gcc
TARGET := solar_system
SRC    := solar_system.c

# ── platform detection ────────────────────────────────────────────────────
UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)
    # Homebrew sets a different prefix on Apple Silicon vs Intel
    BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /usr/local)
    PKG_CONFIG  := PKG_CONFIG_PATH=$(BREW_PREFIX)/lib/pkgconfig pkg-config
else
    PKG_CONFIG  := pkg-config
endif

# ── flags ─────────────────────────────────────────────────────────────────
CFLAGS  := -O2 -std=c17 -Wall -Wextra -Wpedantic \
           -Wno-unused-parameter \
           $(shell $(PKG_CONFIG) --cflags sdl2 SDL2_ttf)

LDFLAGS := $(shell $(PKG_CONFIG) --libs sdl2 SDL2_ttf) -lm

# macOS: also link against the Cocoa framework SDL2 needs
ifeq ($(UNAME), Darwin)
    LDFLAGS += -framework Cocoa
endif

# ── rules ─────────────────────────────────────────────────────────────────
.PHONY: all clean run info

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	$(RM) $(TARGET)

info:
	@echo "Platform : $(UNAME)"
	@echo "Compiler : $(CC)"
	@echo "CFLAGS   : $(CFLAGS)"
	@echo "LDFLAGS  : $(LDFLAGS)"
