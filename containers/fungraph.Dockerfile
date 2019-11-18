ARG IMAGE
ARG FPA_IMAGE
ARG MINIMAP2_IMAGE
ARG MINIGRAPH_IMAGE
ARG GFATOOLS_IMAGE
ARG ODGI_IMAGE
ARG SEQWISH_IMAGE
ARG VG_IMAGE
ARG GRAPHALIGNER_IMAGE
ARG PPG_IMAGE
ARG PYTHON3_IMAGE

FROM "${FPA_IMAGE}" as fpa_builder
FROM "${MINIMAP2_IMAGE}" as minimap2_builder
FROM "${MINIGRAPH_IMAGE}" as minigraph_builder
FROM "${GFATOOLS_IMAGE}" as gfatools_builder
FROM "${ODGI_IMAGE}" as odgi_builder
FROM "${SEQWISH_IMAGE}" as seqwish_builder
FROM "${VG_IMAGE}" as vg_builder
FROM "${GRAPHALIGNER_IMAGE}" as graphaligner_builder
FROM "${PPG_IMAGE}" as ppg_builder
FROM "${PYTHON3_IMAGE}" as python3_builder


FROM "${IMAGE}"

ARG FPA_VERSION
ARG FPA_PREFIX_ARG
ENV FPA_PREFIX="${FPA_PREFIX_ARG}"
LABEL fpa.version="${FPA_VERSION}"

ENV PATH "${FPA_PREFIX}/bin:${PATH}"

COPY --from=fpa_builder "${FPA_PREFIX}" "${FPA_PREFIX}"
COPY --from=fpa_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fpa.txt


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


ARG GFATOOLS_COMMIT
ARG GFATOOLS_PREFIX_ARG
ENV GFATOOLS_PREFIX="${GFATOOLS_PREFIX_ARG}"

LABEL gfatools.version="${GFATOOLS_COMMIT}"

ENV PATH="${GFATOOLS_PREFIX}/bin:${PATH}"

COPY --from=gfatools_builder "${GFATOOLS_PREFIX}" "${GFATOOLS_PREFIX}"
COPY --from=gfatools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/gfatools.txt


ARG MINIGRAPH_COMMIT
ARG MINIGRAPH_REPO
ARG MINIGRAPH_PREFIX_ARG
ENV MINIGRAPH_PREFIX="${MINIGRAPH_PREFIX_ARG}"

LABEL minigraph.version="${MINIGRAPH_COMMIT}"

ENV PATH="${MINIGRAPH_PREFIX}/bin:${PATH}"

COPY --from=minigraph_builder "${MINIGRAPH_PREFIX}" "${MINIGRAPH_PREFIX}"
COPY --from=minigraph_builder "${APT_REQUIREMENTS_FILE}" /build/apt/minigraph.txt


ARG ODGI_COMMIT
ARG ODGI_PREFIX_ARG
ENV ODGI_PREFIX="${ODGI_PREFIX_ARG}"
LABEL odgi.version="${ODGI_VERSION}"

ENV PATH="${ODGI_PREFIX}/bin:${PATH}"

COPY --from=odgi_builder "${ODGI_PREFIX}" "${ODGI_PREFIX}"
COPY --from=odgi_builder "${APT_REQUIREMENTS_FILE}" /build/apt/odgi.txt


ARG SEQWISH_COMMIT
ARG SEQWISH_PREFIX_ARG
ENV SEQWISH_PREFIX="${SEQWISH_PREFIX_ARG}"
LABEL seqwish.version="${SEQWISH_COMMIT}"

ENV PATH "${SEQWISH_PREFIX}/bin:${SEQWISH_COMMIT}/scripts:${PATH}"

COPY --from=seqwish_builder "${SEQWISH_PREFIX}" "${SEQWISH_PREFIX}"
COPY --from=seqwish_builder "${APT_REQUIREMENTS_FILE}" /build/apt/seqwish.txt


ARG VG_VERSION
ARG VG_PREFIX_ARG
ENV VG_PREFIX="${VG_PREFIX_ARG}"

ENV PATH="${VG_PREFIX}/bin:${PATH}"
LABEL vg.version="${VG_VERSION}"

COPY --from=vg_builder "${VG_PREFIX}" "${VG_PREFIX}"
COPY --from=vg_builder "${APT_REQUIREMENTS_FILE}" /build/apt/vg.txt


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

COPY --from=graphaligner_builder "${MUMMER_PREFIX}" "${MUMMER_PREFIX}"


ARG PPG_COMMIT
ARG PPG_PREFIX_ARG
ENV PPG_PREFIX="${PPG_PREFIX_ARG}"

LABEL ppg.version="${PPG_COMMIT}"

ENV PATH "${PPG_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${PPG_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=ppg_builder "${PPG_PREFIX}" "${PPG_PREFIX}"
COPY --from=ppg_builder "${APT_REQUIREMENTS_FILE}" /build/apt/ppg.txt


COPY --from=python3_builder "${APT_REQUIREMENTS_FILE}" /build/apt/python3.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
