PWD=$(shell pwd)
DEBIAN_VERSION=buster-20191014-slim
TAG=fungraph-v0.0.1
MAINTAINER=darcy.ab.jones@gmail.com
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs
MINIMAP2_TAG=v2.17
K8_VERSION=0.2.4
MINIMAP2_REPO=https://github.com/lh3/minimap2.git
K8_URL=https://github.com/attractivechaos/k8/releases/download/v$(K8_VERSION)/k8-$(K8_VERSION).tar.bz2
MINIMAP2_PREFIX_ARG=$(INSTALL_BASE)/minimap2/$(MINIMAP_TAG)
K8_PREFIX_ARG=$(INSTALL_BASE)/k8/$(K8_VERSION)

SEQWISH_COMMIT=62f0055
SEQWISH_REPO=https://github.com/ekg/seqwish.git
SEQWISH_PREFIX_ARG=$(INSTALL_BASE)/seqwish/$(SEQWISH_COMMIT)

FPA_VERSION=0.5.0
FPA_PREFIX_ARG=$(INSTALL_BASE)/fpa/$(FPA_VERSION)

ODGI_COMMIT=73969df
ODGI_REPO=https://github.com/vgteam/odgi.git
ODGI_PREFIX_ARG=$(INSTALL_BASE)/odgi/$(ODGI_COMMIT)

VG_VERSION=v1.19.0
VG_URL=https://github.com/vgteam/vg/releases/download/v1.19.0/vg-v1.19.0.tar.gz
VG_PREFIX_ARG=$(INSTALL_BASE)/vg/$(VG_VERSION)

MINIGRAPH_COMMIT=b1976de
MINIGRAPH_REPO=https://github.com/lh3/minigraph.git
MINIGRAPH_PREFIX_ARG=$(INSTALL_BASE)/minigraph/$(MINIGRAPH_COMMIT)

MUMMER_VERSION=4.0.0beta2
MUMMER_URL=https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
MUMMER_PREFIX_ARG=$(INSTALL_BASE)/mummer/$(MUMMER_VERSION)

GRAPHALIGNER_COMMIT=260df01
GRAPHALIGNER_REPO=https://github.com/maickrau/GraphAligner.git
GRAPHALIGNER_PREFIX_ARG=$(INSTALL_BASE)/graphaligner/$(GRAPHALIGNER_COMMIT)

ZSTR_REPO=https://github.com/mateidavid/zstr.git
BBHASH_REPO=https://github.com/maickrau/BBHash.git
JEMALLOC_URL=https://github.com/jemalloc/jemalloc/releases/download/5.1.0/jemalloc-5.1.0.tar.bz2
SPARSEHASH_URL=https://github.com/sparsehash/sparsehash/archive/sparsehash-2.0.3.tar.gz

GFATOOLS_COMMIT=5f7fe55
GFATOOLS_REPO=https://github.com/lh3/gfatools.git
GFATOOLS_PREFIX_ARG=$(INSTALL_BASE)/gfatools/$(GFATOOLS_COMMIT)

PPG_COMMIT=95cebac
PPG_REPO=https://github.com/darcyabjones/pypafgraph.git
PPG_PREFIX_ARG=$(INSTALL_BASE)/pypafgraph/$(PPG_COMMIT)


GENERAL=base python3 minimap2 seqwish fpa odgi vg minigraph gfatools graphaligner ppg


all: docker/all singularity/all

# Docker targets
docker: docker/panann
docker/all: docker/general

docker/general: $(addprefix docker/,$(GENERAL))

# Singularity targets
singularity: singularity/panann.sif
singularity/all: singularity/general
singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))

# Pull containers for dockerhub
dockerpull: dockerpull/panann
dockerpull/all: dockerpull/general
dockerpull/general: $(addprefix dockerpull/,$(GENERAL))


# Utils
clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean --force
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	singularity cache clean --force
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	@mkdir -p tmp tmp/cache
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --build-arg MAINTAINER="$(MAINTAINER)" \
	  --build-arg PIPELINE_VERSION="$(TAG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/fpa: fpa.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FPA_VERSION="$(FPA_VERSION)" \
	  --build-arg FPA_PREFIX_ARG="$(FPA_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/gfatools: gfatools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GFATOOLS_COMMIT="$(GFATOOLS_COMMIT)" \
	  --build-arg GFATOOLS_REPO="$(GFATOOLS_REPO)" \
	  --build-arg GFATOOLS_PREFIX_ARG="$(GFATOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/minimap2: minimap2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MINIMAP2_TAG="$(MINIMAP2_TAG)" \
	  --build-arg MINIMAP2_REPO="$(MINIMAP2_REPO)" \
	  --build-arg MINIMAP2_PREFIX_ARG="$(MINIMAP2_PREFIX_ARG)" \
	  --build-arg K8_VERSION="$(K8_VERSION)" \
	  --build-arg K8_URL="$(K8_URL)" \
	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/minigraph: minigraph.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MINIGRAPH_COMMIT="$(MINIGRAPH_COMMIT)" \
	  --build-arg MINIGRAPH_REPO="$(MINIGRAPH_REPO)" \
	  --build-arg MINIGRAPH_PREFIX_ARG="$(MINIGRAPH_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/odgi: odgi.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg ODGI_COMMIT="$(ODGI_COMMIT)" \
	  --build-arg ODGI_REPO="$(ODGI_REPO)" \
	  --build-arg ODGI_PREFIX_ARG="$(ODGI_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/seqwish: seqwish.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SEQWISH_COMMIT="$(SEQWISH_COMMIT)" \
	  --build-arg SEQWISH_REPO="$(SEQWISH_REPO)" \
	  --build-arg SEQWISH_PREFIX_ARG="$(SEQWISH_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/vg: vg.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg VG_VERSION="$(VG_VERSION)" \
	  --build-arg VG_URL="$(VG_URL)" \
	  --build-arg VG_PREFIX_ARG="$(VG_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/graphaligner: graphaligner.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GRAPHALIGNER_COMMIT="$(GRAPHALIGNER_COMMIT)" \
	  --build-arg GRAPHALIGNER_REPO="$(GRAPHALIGNER_REPO)" \
	  --build-arg GRAPHALIGNER_PREFIX_ARG="$(GRAPHALIGNER_PREFIX_ARG)" \
	  --build-arg ZSTR_REPO="$(ZSTR_REPO)" \
	  --build-arg BBHASH_REPO="$(BBHASH_REPO)" \
	  --build-arg JEMALLOC_URL="$(JEMALLOC_URL)" \
	  --build-arg SPARSEHASH_URL="$(SPARSEHASH_URL)" \
	  --build-arg MUMMER_VERSION="$(MUMMER_VERSION)" \
	  --build-arg MUMMER_URL="$(MUMMER_URL)" \
	  --build-arg MUMMER_PREFIX_ARG="$(MUMMER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/ppg: ppg.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg PPG_COMMIT="$(PPG_COMMIT)" \
	  --build-arg PPG_REPO="$(PPG_REPO)" \
	  --build-arg PPG_PREFIX_ARG="$(PPG_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/fungraph: docker/fpa docker/minimap2 docker/minigraph docker/odgi
docker/fungraph: docker/seqwish docker/vg docker/graphaligner
docker/fungraph: fungraph.Dockerfile docker/base docker/python3
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FPA_IMAGE="$(call make_image_name,fpa)" \
	  --build-arg FPA_VERSION="$(FPA_VERSION)" \
	  --build-arg FPA_PREFIX_ARG="$(FPA_PREFIX_ARG)" \
	  --build-arg MINIMAP2_IMAGE="$(call make_image_name,minimap2)" \
	  --build-arg MINIMAP2_TAG="$(MINIMAP2_TAG)" \
	  --build-arg MINIMAP2_PREFIX_ARG="$(MINIMAP2_PREFIX_ARG)" \
	  --build-arg K8_VERSION="$(K8_VERSION)" \
	  --build-arg K8_PREFIX_ARG="$(K8_PREFIX_ARG)" \
	  --build-arg MINIGRAPH_IMAGE="$(call make_image_name,minigraph)" \
	  --build-arg MINIGRAPH_COMMIT="$(MINIGRAPH_COMMIT)" \
	  --build-arg MINIGRAPH_PREFIX_ARG="$(MINIGRAPH_PREFIX_ARG)" \
	  --build-arg GFATOOLS_IMAGE="$(call make_image_name,gfatools)" \
	  --build-arg GFATOOLS_COMMIT="$(GFATOOLS_COMMIT)" \
	  --build-arg GFATOOLS_PREFIX_ARG="$(GFATOOLS_PREFIX_ARG)" \
	  --build-arg ODGI_IMAGE="$(call make_image_name,odgi)" \
	  --build-arg ODGI_COMMIT="$(ODGI_COMMIT)" \
	  --build-arg ODGI_PREFIX_ARG="$(ODGI_PREFIX_ARG)" \
	  --build-arg SEQWISH_IMAGE="$(call make_image_name,seqwish)" \
	  --build-arg SEQWISH_COMMIT="$(SEQWISH_COMMIT)" \
	  --build-arg SEQWISH_PREFIX_ARG="$(SEQWISH_PREFIX_ARG)" \
	  --build-arg VG_IMAGE="$(call make_image_name,vg)" \
	  --build-arg VG_VERSION="$(VG_VERSION)" \
	  --build-arg VG_PREFIX_ARG="$(VG_PREFIX_ARG)" \
	  --build-arg GRAPHALIGNER_IMAGE="$(call make_image_name,graphaligner)" \
	  --build-arg GRAPHALIGNER_COMMIT="$(GRAPHALIGNER_COMMIT)" \
	  --build-arg GRAPHALIGNER_PREFIX_ARG="$(GRAPHALIGNER_PREFIX_ARG)" \
	  --build-arg MUMMER_VERSION="$(MUMMER_VERSION)" \
	  --build-arg MUMMER_PREFIX_ARG="$(MUMMER_PREFIX_ARG)" \
	  --build-arg PPG_IMAGE="$(call make_image_name,ppg)" \
	  --build-arg PPG_COMMIT="$(PPG_COMMIT)" \
	  --build-arg PPG_PREFIX_ARG="$(PPG_PREFIX_ARG)" \
	  --build-arg PYTHON3_IMAGE="$(call make_image_name,python3)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@
