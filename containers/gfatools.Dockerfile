ARG IMAGE

FROM "${IMAGE}" as gfatools_builder

ARG GFATOOLS_COMMIT
ARG GFATOOLS_REPO
ARG GFATOOLS_PREFIX_ARG
ENV GFATOOLS_PREFIX="${GFATOOLS_PREFIX_ARG}"


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
  && git clone "${GFATOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "${GFATOOLS_COMMIT}" \
  && make \
  && mkdir -p "${GFATOOLS_PREFIX}/bin" \
  && cp gfatools "${GFATOOLS_PREFIX}/bin" \
  && add_runtime_dep zlib1g


FROM "${IMAGE}"

ARG GFATOOLS_COMMIT
ARG GFATOOLS_REPO
ARG GFATOOLS_PREFIX_ARG
ENV GFATOOLS_PREFIX="${GFATOOLS_PREFIX_ARG}"

LABEL gfatools.version="${GFATOOLS_COMMIT}"

ENV PATH="${GFATOOLS_PREFIX}/bin:${PATH}"

COPY --from=gfatools_builder "${GFATOOLS_PREFIX}" "${GFATOOLS_PREFIX}"
COPY --from=gfatools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gfatools.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
