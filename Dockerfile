FROM ubuntu:16.04

RUN apt-get update \
    && apt-get install -yq --no-install-recommends \
    build-essential ca-certificates libgmp3-dev git software-properties-common python-software-properties \
    wget unzip zlib1g-dev python python-dev python-setuptools parallel samtools maven libgmp-dev gettext automake autopoint libtool libtre-dev \
    && easy_install pip && pip install numpy pandas matplotlib gmpy pairwise pysam biopython tre

# install java
RUN \
    echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-get update && \
    apt-get install -y oracle-java8-installer

# install STAR and art_illumina
RUN cd /opt \
    && wget --quiet -O star.zip https://github.com/alexdobin/STAR/archive/2.5.3a.zip \
 	  && unzip star.zip \
 	  && make -C STAR-2.5.3a/source \
 	  && wget --quiet -O art_illumina.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz \
    && tar -xzvf art_illumina.tgz \
    && wget --quiet https://github.com/repseqio/repseqio/releases/download/v1.2.8/repseqio-1.2.8.zip \
    && unzip repseqio-1.2.8.zip \
    && wget --quiet https://github.com/milaboratory/mitools/releases/download/v1.5/mitools-1.5.zip \
    && unzip mitools-1.5.zip \
    && wget --quiet https://github.com/milaboratory/mixcr/releases/download/v2.1.3/mixcr-2.1.3.zip \
    && unzip mixcr-2.1.3.zip \
    && wget --quiet http://www.nature.com/ng/journal/v49/n4/extref/ng.3820-s2.zip \
    && unzip ng.3820-s2.zip \
    && mv SupplementarySoftware TRUST \
    && rm -rf TRUST/.git \
    && rm -r -f __MACOSX 2>/dev/null || true \
    && rm *.zip *.tgz

ENV PATH="/opt/mixcr-2.1.3:/opt/mitools-1.5:/opt/repseqio-1.2.8:/opt/scripts:/opt/art_bin_MountRainier:/opt/STAR-2.5.3a/source:${PATH}"

ADD falseExtensionsStat.py falseOverlapsStat.py getFalseExtensions.py getFalseOverlaps.py plotMiXCRvsTRUST.py run-comparison.sh run-false-positives.sh /opt/scripts/ 
WORKDIR /work

# ENTRYPOINT /opt/scripts/run-comparison.sh