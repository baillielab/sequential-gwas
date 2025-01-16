# sequential-gwas
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://baillielab.github.io/sequential-gwas/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://baillielab.github.io/sequential-gwas/dev/)
[![Build Status](https://github.com/baillielab/sequential-gwas/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/baillielab/sequential-gwas/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/baillielab/sequential-gwas/branch/main/graph/badge.svg)](https://codecov.io/gh/baillielab/sequential-gwas)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Build

```bash
docker build -t sequential-gwas -f docker/Dockerfile .
```

If running on MacOS with arm platform, add: `--platform linux/amd64`

## Usage

### Regenie (via Docker)

```bash
docker run -it --rm sequential-gwas /root/miniforge3/bin/mamba run -n regenie_env regenie --help
```

If running on MacOS with arm platform, add: `--platform linux/amd64`