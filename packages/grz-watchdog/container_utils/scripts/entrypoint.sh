#!/bin/bash
set -euo pipefail

pixi run --manifest-path /grz-watchdog/pixi.toml -- "$@"
