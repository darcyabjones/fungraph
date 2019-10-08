ARG IMAGE

FROM "${IMAGE}" as odgi_builder

ARG ODGI_COMMIT
ARG ODGI_REPO
ARG ODGI_PREFIX_ARG
ENV ODGI_PREFIX="${ODGI_PREFIX_ARG}"

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
  && git clone "${ODGI_REPO}" . \
  && cmake -DBUILD_STATIC=1 -H. -Bbuild \
  && cmake --build build -- -j 3 \
  && mkdir -p "${ODGI_PREFIX}/bin" \
  && cp bin/odgi "${ODGI_PREFIX}/bin" \
  && add_runtime_dep libgomp1 zlib1g


FROM "${IMAGE}"

ARG ODGI_COMMIT
ARG ODGI_PREFIX_ARG
ENV ODGI_PREFIX="${ODGI_PREFIX_ARG}"
LABEL odgi.version="${ODGI_VERSION}"

ENV PATH="${ODGI_PREFIX}/bin:${PATH}"

COPY --from=odgi_builder "${ODGI_PREFIX}" "${ODGI_PREFIX}"
COPY --from=odgi_builder "${APT_REQUIREMENTS_FILE}" /build/apt/odgi.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
