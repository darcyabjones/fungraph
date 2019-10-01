ARG IMAGE

FROM "${IMAGE}" as fpa_builder

ARG FPA_VERSION
ARG FPA_PREFIX_ARG
ENV FPA_PREFIX="${FPA_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       ca-certificates \
       cargo \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && cargo install --version "${FPA_VERSION}" --root "${FPA_PREFIX}" fpa_lr \
  && add_runtime_dep libbz2-1.0 lzma zlib1g


FROM "${IMAGE}"

ARG FPA_VERSION
ARG FPA_PREFIX_ARG
ENV FPA_PREFIX="${FPA_PREFIX_ARG}"
LABEL fpa.version="${FPA_VERSION}"

ENV PATH "${FPA_PREFIX}/bin:${PATH}"

COPY --from=fpa_builder "${FPA_PREFIX}" "${FPA_PREFIX}"
COPY --from=fpa_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fpa.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
