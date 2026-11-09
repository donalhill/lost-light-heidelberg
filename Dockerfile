# The one environment, the same on the laptop (under Rosetta) and on AWS. Ubuntu 24.04 and its Python 3.12.3, as the
# dry runs' instances, with the packages heidelberg_dry3's instance installed (requirements.lock, its job_requirements.txt).
FROM --platform=linux/amd64 ubuntu:24.04@sha256:fe81560990fae86dc4743202b2cb4d3df98d0bcf949173cecef15ce1127ef852
ENV DEBIAN_FRONTEND=noninteractive PYTHONDONTWRITEBYTECODE=1
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv libcfitsio-bin ca-certificates curl git \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.lock /tmp/requirements.lock
RUN python3 -m venv /opt/v && /opt/v/bin/pip install --no-cache-dir -r /tmp/requirements.lock
RUN python3 -m venv /opt/sm && /opt/sm/bin/pip install --no-cache-dir snakemake==9.16.3 awscli==1.46.2 \
    && ln -s /opt/sm/bin/snakemake /opt/sm/bin/aws /usr/local/bin/
COPY tools/tlock/linux_amd64/tle /usr/local/bin/tle
# numpy's linear algebra on one thread per process: results that do not depend on the machine's size (its sums are
# otherwise split across every core), and no oversubscription when many pairs or toys run as separate processes
ENV PATH=/opt/v/bin:$PATH OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKDIR /repo
