ARG IMAGE


FROM "${IMAGE}" as mummer_builder

ARG MUMMER_VERSION
ARG MUMMER_URL
ARG MUMMER_PREFIX_ARG
ENV MUMMER_PREFIX="${MUMMER_PREFIX_ARG}"

ENV PATH="${MUMMER_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${MUMMER_PREFIX}/lib:${PATH}"
ENV CPATH="${MUMMER_PREFIX}/include:${CPATH}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O mummer.tar.gz "${MUMMER_URL}" \
  && tar -zxf mummer.tar.gz \
  && cd mummer-*/ \
  && ./configure --prefix="${MUMMER_PREFIX}" \
  && make \
  && make install \
  && add_runtime_dep perl sed gawk fig2dev gnuplot-nox xfig


FROM "${IMAGE}" as graphaligner_builder

ARG GRAPHALIGNER_COMMIT
ARG GRAPHALIGNER_REPO
ARG GRAPHALIGNER_PREFIX_ARG
ENV GRAPHALIGNER_PREFIX="${GRAPHALIGNER_PREFIX_ARG}"

ARG ZSTR_REPO
ARG BBHASH_REPO
ARG JEMALLOC_URL
ARG SPARSEHASH_URL

ARG MUMMER_VERSION
ARG MUMMER_PREFIX_ARG
ENV MUMMER_PREFIX="${MUMMER_PREFIX_ARG}"

ENV PATH="${MUMMER_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${MUMMER_PREFIX}/lib:${PATH}"
ENV PKG_CONFIG_PATH="${MUMMER_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}"
ENV CPATH="${MUMMER_PREFIX}/include:${CPATH}"

COPY --from=mummer_builder "${MUMMER_PREFIX}" "${MUMMER_PREFIX}"
COPY --from=mummer_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mummer.txt

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       libboost-dev \
       libboost-program-options-dev \
       libboost-serialization-dev \
       libjemalloc-dev \
       libprotobuf-dev \
       libsdsl-dev \
       pkg-config \
       protobuf-compiler \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${GRAPHALIGNER_REPO}" . \
  && rm -rf -- BBHash zstr \
  && git clone "${ZSTR_REPO}" \
  && git clone "${BBHASH_REPO}" \
  && git submodule update --init --recursive \
  && wget "${JEMALLOC_URL}" \
  && tar -xjf jemalloc*.tar.bz2 \
  && cd jemalloc*/ \
  && ./configure && make && make install \
  && cd .. \
  && wget "${SPARSEHASH_URL}" \
  && tar -zxf sparsehash*.tar.gz \
  && cd sparsehash*/ \
  && ./configure && make && make install \
  && cd .. \
  && make bin/GraphAligner \
  && mkdir -p "${GRAPHALIGNER_PREFIX}/bin" \
  && cp bin/GraphAligner "${GRAPHALIGNER_PREFIX}/bin" \
  && add_runtime_dep \
       libprotobuf17 \
       libgomp1 \
       libprotoc17 \
       libsdsl3 \
       libjemalloc2 \
       libboost-program-options1.67.0


FROM "${IMAGE}"

ARG GRAPHALIGNER_COMMIT
ARG GRAPHALIGNER_PREFIX_ARG
ENV GRAPHALIGNER_PREFIX="${GRAPHALIGNER_PREFIX_ARG}"
LABEL graphaligner.version="${GRAPHALIGNER_COMMIT}"

ENV PATH="${GRAPHALIGNER_PREFIX}/bin:${PATH}"

COPY --from=graphaligner_builder "${GRAPHALIGNER_PREFIX}" "${GRAPHALIGNER_PREFIX}"
COPY --from=graphaligner_builder "${APT_REQUIREMENTS_FILE}" /build/apt/graphaligner_builder.txt

ARG MUMMER_VERSION
ARG MUMMER_PREFIX_ARG
ENV MUMMER_PREFIX="${MUMMER_PREFIX_ARG}"

LABEL mummer.version="${MUMMER_VERSION}"

ENV PATH="${MUMMER_PREFIX}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${MUMMER_PREFIX}/lib:${PATH}"
ENV PKG_CONFIG_PATH="${MUMMER_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}"
ENV CPATH="${MUMMER_PREFIX}/include:${CPATH}"

COPY --from=mummer_builder "${MUMMER_PREFIX}" "${MUMMER_PREFIX}"
COPY --from=mummer_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mummer.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
