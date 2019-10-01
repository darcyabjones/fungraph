ARG IMAGE

FROM "${IMAGE}" as builder

ARG SEQTK_COMMIT
ARG SEQTK_REPO="https://github.com/lh3/seqtk.git"
ARG SEQTK_PREFIX_ARG="/opt/seqtk/${SEQTK_COMMIT}"
ENV SEQTK_PREFIX="${SEQTK_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SEQTK_REPO}" . \
  && git fetch --tags \
  && git checkout "${SEQTK_COMMIT}" \
  && make \
  && mkdir -p "${SEQTK_PREFIX}/bin" \
  && make BINDIR="${SEQTK_PREFIX}/bin" install \
  && add_runtime_dep zlib1g


FROM "${IMAGE}"

ARG SEQTK_COMMIT
ARG SEQTK_PREFIX_ARG="/opt/seqtk/${SEQTK_COMMIT}"
ENV SEQTK_PREFIX="${SEQTK_PREFIX_ARG}"
LABEL seqtk.version="${SEQTK_COMMIT}"

ENV PATH "${SEQTK_PREFIX}/bin:${PATH}"

COPY --from=builder "${SEQTK_PREFIX}" "${SEQTK_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/seqtk.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
