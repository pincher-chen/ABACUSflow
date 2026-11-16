# syntax=docker/dockerfile:1
#
# ABACUSflow image -- based on the official ABACUS 3.11 GNU (OpenMPI) image,
# with the ABACUSflow CLI (workflow / input generation) layered on top.
#
#   base : ghcr.io/deepmodeling/abacus-gnu:3.11.0-beta4
#          Ubuntu 22.04.5, Python 3.10.12, ABACUS v3.11.0-beta.4
#          (/usr/local/bin/abacus -> abacus_max_para), OpenMPI 4.1.2.
#          Environment already sets OMPI_ALLOW_RUN_AS_ROOT(_CONFIRM) so the
#          container can run MPI as root without --allow-run-as-root.
#
# Difference from the Intel build: no oneAPI, no MKL; the image is ~2.4 GB
# instead of ~9.9 GB. The base has no python/pip aliases and no
# ca-certificates, so we create the `python` symlink ourselves and let uv use
# its bundled root store for PyPI over TLS.
#
# Smoke test:
#   python -m abacusflow.abacus --help
#   python /opt/abacusflow/abacus.py --help
#
# Build (from the directory containing this Dockerfile and ABACUSflow/):
#   podman build -t abacusflow:3.11-gnu -f Dockerfile.3.11-gnu .
FROM ghcr.io/deepmodeling/abacus-gnu:3.11.0-beta4

LABEL org.opencontainers.image.title="abacusflow-gnu" \
      org.opencontainers.image.description="ABACUSflow workflow generator on top of the ABACUS 3.11 GNU/OpenMPI runtime (abacus-gnu:3.11.0-beta4)" \
      org.opencontainers.image.source="https://github.com/pincher-chen/ABACUSflow.git" \
      org.opencontainers.image.base.name="ghcr.io/deepmodeling/abacus-gnu:3.11.0-beta4"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    PATH=/opt/abacusflow/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    # `python -m abacusflow.abacus` must work from any directory
    PYTHONPATH=/opt/abacusflow

# --- Python tooling --------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

RUN set -eux; \
    # The base image ships python3 only; the CLI and the generated job scripts
    # invoke `python`.
    ln -sf /usr/bin/python3 /usr/local/bin/python; \
    python -V

# --- Python dependencies ---------------------------------------------------
COPY requirements.txt /tmp/requirements.txt
RUN UV_PYTHON_DOWNLOADS=never uv pip install --system --python "$(command -v python)" --no-cache -r /tmp/requirements.txt \
 && rm -f /tmp/requirements.txt \
 && python -c "import ase, numpy, pymatgen, click, yaml; print('deps ok')"

# --- ABACUSflow source -----------------------------------------------------
COPY ABACUSflow/ /opt/abacusflow/

RUN set -eux; \
    find /opt/abacusflow -name '__pycache__' -type d -prune -exec rm -rf {} +; \
    # Stale local backup configs must never ship: they carry the original
    # author's cluster paths and could act as a fallback configuration.
    find /opt/abacusflow/config -maxdepth 1 -name 'condor.ini.bak*' -delete; \
    rm -f /opt/abacusflow/config/condor.ini.example; \
    rm -rf /opt/abacusflow/.git; \
    \
    # Upstream is a flat tree (abacus.py, abacus/, config/). Wrap it in an
    # `abacusflow` package while keeping the upstream layout inside, so
    # `from abacus.xxx import ...` and the `from config import ...` fallback
    # keep resolving. `abacus/__main__.py` makes `python -m abacusflow.abacus`
    # work (in CPython the directory `abacus/` shadows sibling `abacus.py`).
    mkdir -p /opt/abacusflow/abacusflow; \
    cp -a /opt/abacusflow/abacus.py /opt/abacusflow/abacus /opt/abacusflow/config /opt/abacusflow/abacusflow/; \
    printf '"""ABACUSflow: ABACUS workflow generation toolkit."""\n\n__version__ = "0.1.0"\n' \
        > /opt/abacusflow/abacusflow/__init__.py; \
    printf '%s\n' \
        '"""Entry point for `python -m abacusflow.abacus`.' \
        '' \
        'ABACUSflow is a flat script tree (abacus.py + abacus/ + config/), so in' \
        'CPython the `abacus` subpackage shadows the sibling abacus.py module for' \
        'this dotted name. Delegate to the CLI defined in abacusflow/abacus.py.' \
        '"""' \
        '' \
        'import os as _os' \
        '' \
        '_CLI_SCRIPT = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "abacus.py")' \
        '' \
        'if __name__ == "__main__":' \
        '    import runpy' \
        '' \
        '    runpy.run_path(_CLI_SCRIPT, run_name="__main__")' \
        > /opt/abacusflow/abacusflow/abacus/__main__.py; \
    \
    mkdir -p /opt/abacusflow/bin; \
    printf '#!/bin/sh\nexec python -m abacusflow.abacus "$@"\n' > /opt/abacusflow/bin/abacusflow; \
    printf '#!/bin/sh\nexec python /opt/abacusflow/abacus.py "$@"\n' > /opt/abacusflow/bin/abacus.py; \
    chmod 0755 /opt/abacusflow/bin/abacusflow /opt/abacusflow/bin/abacus.py; \
    test -f /opt/abacusflow/abacusflow/abacus.py;

# --- Smoke test at build time ---------------------------------------------
RUN cd / && python -m abacusflow.abacus --help > /dev/null \
 && cd /opt/abacusflow && python abacus.py --help > /dev/null \
 && abacus --version > /dev/null \
 && rm -rf /opt/abacusflow/abacusflow/__pycache__ /opt/abacusflow/abacusflow/abacus/__pycache__

# --- Non-root user and working dir ----------------------------------------
RUN useradd --create-home --uid 1000 abacus \
 && mkdir -p /work && chown abacus:abacus /work
USER abacus

WORKDIR /opt/abacusflow
CMD ["python", "abacus.py"]
