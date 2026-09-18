#!/bin/zsh
cd "$(dirname "$0")" || exit 1
exec .venv/bin/jupyter lab notebooks/01_static_mzm.ipynb
