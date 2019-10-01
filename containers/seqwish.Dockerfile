ARG IMAGE

FROM "${IMAGE}" as seqwish_builder

ARG SEQWISH_COMMIT
ARG SEQWISH_REPO
ARG SEQWISH_PREFIX_ARG
ENV SEQWISH_PREFIX="${SEQWISH_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       git \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${SEQWISH_REPO}" . \
  && cmake -H. -Bbuild \
  && cmake --build build -- -j3 \
  && mkdir -p "${SEQWISH_PREFIX}" \
  && cp -r "bin" "${SEQWISH_PREFIX}" \
  && cp -r "scripts" "${SEQWISH_PREFIX}" \
  && add_runtime_dep zlib1g libgomp1 libatomic1


FROM "${IMAGE}"

ARG SEQWISH_COMMIT
ARG SEQWISH_PREFIX_ARG
ENV SEQWISH_PREFIX="${SEQWISH_PREFIX_ARG}"
LABEL seqwish.version="${SEQWISH_COMMIT}"

ENV PATH "${SEQWISH_PREFIX}/bin:${SEQWISH_COMMIT}/scripts:${PATH}"

COPY --from=seqwish_builder "${SEQWISH_PREFIX}" "${SEQWISH_PREFIX}"
COPY --from=seqwish_builder "${APT_REQUIREMENTS_FILE}" /build/apt/seqwish.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
