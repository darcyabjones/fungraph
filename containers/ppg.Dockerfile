ARG IMAGE
FROM "${IMAGE}" as ppg_builder

ARG PPG_COMMIT
ARG PPG_REPO
ARG PPG_PREFIX_ARG
ENV PPG_PREFIX="${PPG_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       python3 \
       python3-pip \
       python3-setuptools \
       python3-wheel \
       python3-intervaltree \
       python3-biopython \
       python3-networkx \
       python3-scipy \
       python3-matplotlib \
       python3-tk \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${PPG_REPO}" . \
  && git fetch --tags \
  && git checkout "${PPG_COMMIT}" \
  && pip3 install --prefix="${PPG_PREFIX}" . \
  && add_runtime_dep \
       python3 \
       python3-intervaltree \
       python3-biopython \
       python3-networkx \
       python3-scipy \
       python3-matplotlib \
       python3-tk


FROM "${IMAGE}"

ARG PPG_COMMIT
ARG PPG_PREFIX_ARG
ENV PPG_PREFIX="${PPG_PREFIX_ARG}"

LABEL ppg.version="${PPG_COMMIT}"

ENV PATH "${PPG_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${PPG_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=ppg_builder "${PPG_PREFIX}" "${PPG_PREFIX}"
COPY --from=ppg_builder "${APT_REQUIREMENTS_FILE}" /build/apt/ppg.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
