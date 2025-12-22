#!/usr/bin/env bash
set -euo pipefail
CCX="${CCX:-ccx}"
PY="${PY:-python3}"
cd "$(dirname "$0")"

echo "[beam_patch] Test 1: cantilever"
$CCX -i cantilever_b32r_unitload
$PY check_dat.py cantilever_b32r_unitload.dat --node 1 --comp U1 --expected 0.02250649 --rtol 0.10

echo "[beam_patch] Test 2: simply supported"
$CCX -i simply_supported_b32r_midload
$PY check_dat.py simply_supported_b32r_midload.dat --node 2 --comp U1 --expected 0.00113403 --rtol 0.10

echo "[beam_patch] ALL OK"
