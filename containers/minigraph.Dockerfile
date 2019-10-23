ARG IMAGE

FROM "${IMAGE}" as minigraph_builder

ARG MINIGRAPH_COMMIT
ARG MINIGRAPH_REPO
ARG MINIGRAPH_PREFIX_ARG
ENV MINIGRAPH_PREFIX="${MINIGRAPH_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       curl \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${MINIGRAPH_REPO}" . \
  && git fetch --tags \
  && git checkout "${MINIGRAPH_COMMIT}" \
  && make \
  && mkdir -p "${MINIGRAPH_PREFIX}/bin" \
  && cp minigraph "${MINIGRAPH_PREFIX}/bin" \
  && add_runtime_dep zlib1g


FROM "${IMAGE}"

ARG MINIGRAPH_COMMIT
ARG MINIGRAPH_REPO
ARG MINIGRAPH_PREFIX_ARG
ENV MINIGRAPH_PREFIX="${MINIGRAPH_PREFIX_ARG}"

LABEL minigraph.version="${MINIGRAPH_COMMIT}"

ENV PATH="${MINIGRAPH_PREFIX}/bin:${PATH}"

COPY --from=minigraph_builder "${MINIGRAPH_PREFIX}" "${MINIGRAPH_PREFIX}"
COPY --from=minigraph_builder "${APT_REQUIREMENTS_FILE}" /build/apt/minigraph.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
