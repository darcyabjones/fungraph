ARG IMAGE

FROM "${IMAGE}" as minimap2_builder

ARG MINIMAP2_TAG
ARG MINIMAP2_REPO="https://github.com/lh3/minimap2.git"
ARG MINIMAP2_PREFIX_ARG="/opt/minimap2/${MINIMAP2_TAG}"
ENV MINIMAP2_PREFIX="${MINIMAP2_PREFIX_ARG}"

ARG K8_VERSION
ARG K8_URL="https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2"
ARG K8_PREFIX_ARG="/opt/k8/${K8_VERSION}"
ENV K8_PREFIX="${K8_PREFIX_ARG}"

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
  && git clone "${MINIMAP2_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${MINIMAP2_TAG}" \
  && make extra \
  && mkdir -p "${MINIMAP2_PREFIX}/bin" \
  && cp minimap2 sdust minimap2-lite "${MINIMAP2_PREFIX}/bin" \
  && cp misc/paftools.js "${MINIMAP2_PREFIX}/bin" \
  && curl -L "${K8_URL}" \
   | tar jxf - \
  && mkdir -p "${K8_PREFIX}/bin" \
  && mv "k8-${K8_VERSION}/k8-$(uname -s)" "${K8_PREFIX}/bin/k8" \
  && add_runtime_dep zlib1g

# CA cert stuff sometime required for git clone https


FROM "${IMAGE}"

ARG MINIMAP2_TAG
ARG MINIMAP2_PREFIX_ARG="/opt/minimap2/${MINIMAP2_TAG}"
ENV MINIMAP2_PREFIX="${MINIMAP2_PREFIX_ARG}"
LABEL minimap2.version="${MINIMAP2_TAG}"

ARG K8_VERSION
ARG K8_PREFIX_ARG="/opt/k8/${K8_VERSION}"
ENV K8_PREFIX="${K8_PREFIX_ARG}"
LABEL k8.version="${K8_VERSION}"

ENV PATH "${MINIMAP2_PREFIX}/bin:${K8_PREFIX}/bin:${PATH}"

COPY --from=minimap2_builder "${MINIMAP2_PREFIX}" "${MINIMAP2_PREFIX}"
COPY --from=minimap2_builder "${K8_PREFIX}" "${K8_PREFIX}"
COPY --from=minimap2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/minimap2.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
