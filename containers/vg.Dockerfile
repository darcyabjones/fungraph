ARG IMAGE

FROM "${IMAGE}" as vg_builder

ARG VG_VERSION
ARG VG_URL
ARG VG_PREFIX_ARG
ENV VG_PREFIX="${VG_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       automake \
       bc \
       bison \
       build-essential \
       ca-certificates \
       cmake \
       curl \
       flex \
       git \
       gawk \
       jq \
       libbz2-dev \
       libcairo2-dev \
       libffi-dev \
       libjansson-dev \
       liblz4-dev \
       liblzma-dev \
       libncurses-dev \
       libprotobuf-dev \
       libprotoc-dev \
       librdf0-dev \
       libtool \
       lzma-dev \
       pkg-config \
       protobuf-compiler \
       redland-utils \
       unzip \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O vg.tar.gz "${VG_URL}" \
  && tar -zxf vg.tar.gz \
  && cd vg-*/ \
  && make \
  && mkdir -p "${VG_PREFIX}" \
  && cp -r bin "${VG_PREFIX}" \
  && add_runtime_dep \
       libatomic1 \
       libboost-program-options1.67.0 \
       libbz2-1.0 \
       libcairo2 \
       libcairo-gobject2 \
       libffi6 \
       libgomp1 \
       libjansson4 \
       libjemalloc2 \
       liblz4-1 \
       liblzma5 \
       libncurses6 \
       libprotobuf17 \
       libprotoc17 \
       librdf0 \
       zlib1g


FROM "${IMAGE}"

ARG VG_VERSION
ARG VG_PREFIX_ARG
ENV VG_PREFIX="${VG_PREFIX_ARG}"

ENV PATH="${VG_PREFIX}/bin:${PATH}"
LABEL vg.version="${VG_VERSION}"

COPY --from=vg_builder "${VG_PREFIX}" "${VG_PREFIX}"
COPY --from=vg_builder "${APT_REQUIREMENTS_FILE}" /build/apt/vg.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
