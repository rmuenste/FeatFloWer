# Container image for building and running FeatFloWer
FROM ubuntu:24.04 AS base

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      g++ \
      gcc \
      cmake \
      git \
      ninja-build \
      openmpi-bin \
      libopenmpi-dev \
      liblapack-dev \
      libblas-dev \
      libmetis-dev \
      libboost-dev \
      libboost-program-options-dev \
      ca-certificates \
      python3 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Copy source tree into the image; expect submodules alongside the checkout.
COPY . /workspace

RUN git submodule update --init --recursive

# Configure and build once during image creation so downstream users can run binaries directly.
ARG CMAKE_BUILD_TYPE=Release
ARG BUILD_APPLICATIONS=ON
ARG Q2P1_BUILD_ID=generic-linux-gcc-release
RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
        -DBUILD_APPLICATIONS=${BUILD_APPLICATIONS} \
        -DQ2P1_BUILD_ID=${Q2P1_BUILD_ID} && \
    cmake --build build -- -j"$(nproc)"

# Provide a non-root user to simplify volume mounts on HPC/cloud systems.
ARG USERNAME=featflower
ARG USER_UID=1000
ARG USER_GID=1000
RUN set -eux; \
    existing_group="$(getent group "${USER_GID}" | cut -d: -f1 || true)"; \
    if [ -z "${existing_group}" ]; then \
      groupadd --gid "${USER_GID}" "${USERNAME}"; \
    elif [ "${existing_group}" != "${USERNAME}" ]; then \
      groupmod -n "${USERNAME}" "${existing_group}"; \
    fi; \
    existing_user="$(getent passwd "${USER_UID}" | cut -d: -f1 || true)"; \
    if [ -z "${existing_user}" ]; then \
      if ! id -u "${USERNAME}" >/dev/null 2>&1; then \
        useradd --uid "${USER_UID}" --gid "${USER_GID}" --shell /bin/bash --create-home "${USERNAME}"; \
      else \
        usermod --uid "${USER_UID}" --gid "${USER_GID}" --shell /bin/bash "${USERNAME}"; \
      fi; \
    elif [ "${existing_user}" != "${USERNAME}" ]; then \
      usermod --login "${USERNAME}" --home "/home/${USERNAME}" --move-home "${existing_user}"; \
      usermod --uid "${USER_UID}" --gid "${USER_GID}" --shell /bin/bash "${USERNAME}"; \
    else \
      usermod --uid "${USER_UID}" --gid "${USER_GID}" --shell /bin/bash "${USERNAME}"; \
    fi

USER ${USERNAME}
WORKDIR /workspace

# Default to a shell so users can run applications, tests, or MPI jobs via docker run.
CMD ["/bin/bash"]
