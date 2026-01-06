# ItsNotK1C

**ItsNotK1C** is a heavily rebuilt Klipper-based 3D printer platform that
originated from a Creality K1C chassis but no longer follows any stock K1/K1C
architectural, electrical, or mechanical assumptions.

At this point, the original K1C identity is purely historical.  
This machine should be treated as a custom Klipper integration and research
platform rather than a modified consumer printer.

---

## Project Links

- **Klipper (base firmware)** https://www.klipper3d.org/
- **Kalico (community-maintained fork firmware)** https://github.com/KalicoCrew/kalico
- **SimpleAF project (supporting context)** https://pellcorp.github.io/creality-wiki/
- **BoxTurtle MMU (reference design)** https://www.armoredturtle.xyz/docs/boxturtle/index.html


## High-level Overview

- Klipper-based, multi-MCU architecture
- External host offloading using Raspberry Pi 5
- Fully rebuilt motion system
- AC-powered heated bed with tooling plate
- Custom toolhead, hotend, and extruder
- Custom MMU derived from BoxTurtle
- Multiple probing systems

The system exists to overcome host SoC limitations, timing instability
(e.g. *timer too close*), and architectural constraints inherent to the original
K1 platform.

---

## Software Stack

- **Base firmware**: Mainstream Klipper
- **Codebase strategy**:
  - Derived from upstream Klipper
  - Tracks mainstream Klipper development for long-term rebaseability
  - Avoids vendor-specific assumptions tied to K1/K1C firmware

- **Architecture focus**:
  - External host offloading
  - Multi-MCU coordination
  - Deterministic scheduling and timing stability
  - Support for non-standard hardware topology

- **Project support**:
  - This work is supported by the **SimpleAF project**
  - SimpleAF provides the broader experimentation context and long-term
    architectural direction for this platform

The goal of the software stack is architectural correctness and system stability,
not feature accumulation.

---

## Host & Control Architecture

- **Primary control firmware**: Klipper (customized)
- **Host topology**:
  - Original K1 SoC retained only for minimal platform functions
  - Raspberry Pi 5 used as an **external Klipper host**
  - Host-to-host communication bridged via `socat`

- **Motivation**:
  - Eliminate scheduling instability caused by the weak stock SoC
  - Provide sufficient compute headroom for multi-MCU operation,
    MMU logic, and advanced probing

---

## MCU Topology

- Multi-MCU configuration
- Independent MCUs assigned to:
  - Main motion control
  - Toolhead
  - Bed and power-related peripherals
  - Auxiliary subsystems (MMU, probes, etc.)

The MCU layout is designed to scale beyond the assumptions of the original
K1/K1C electronics.

---

## Motion System

- All stock motion components removed
- Linear rails on X, Y, and Z axes
- Triple independent Z axes
- Oldham couplers
- All stepper motors replaced with **LDO steppers**

The resulting motion system shares no functional similarity with the stock K1/K1C
implementation.

---

## Bed System

- **AC-powered heated bed (220V)**
- Solid-state relay (SSR) controlled
- **MIC6 aluminum tooling plate**
  - Thickness: 6 mm
- Complete replacement of the original 24V DC bed system

---

## Toolhead & Extrusion

- Fully custom toolhead
- Non-stock hotend
- Non-stock extruder
- Designed to operate correctly within:
  - Multi-MCU architecture
  - External host offload
  - Custom MMU integration
  - Advanced probing workflows

No original K1/K1C toolhead assumptions apply.

---

## MMU (Multi-Material Unit)

- Custom MMU implementation
- **Derived from BoxTurtle**, but:
  - Mechanically modified
  - Electrically modified
  - Firmware and configuration heavily adapted

This MMU is **not compatible** with stock BoxTurtle configurations.

---

## Probing & Calibration

- Multiple probing systems deployed
- Combination of:
  - Contact-based probing
  - Non-contact / eddy-style probing

Used for:
- Bed leveling
- Z calibration
- Redundancy and experimentation

---

## Naming Rationale

**ItsNotK1C** is an intentional name.

This machine:
- Is not electrically a K1C
- Is not mechanically a K1C
- Is not architecturally a K1C
- Does not conform to stock K1/K1C firmware, timing, or control assumptions

The name reflects current reality, not original branding.

---

## Intended Use

- Klipper experimentation
- Multi-MCU and host offloading research
- MMU development and testing
- Timing, scheduling, and reliability validation
- Integration of non-standard hardware combinations

This repository is **not intended** to be a drop-in solution for stock K1/K1C users.

---

## Status

- Actively modified
- Architecture-first, not vendor-first
- No guarantees of compatibility with any Creality firmware or hardware

---

## Disclaimer

This project is based on Klipper and integrates ideas and code from multiple
downstream forks and experiments.

It is provided strictly for research and development purposes.
