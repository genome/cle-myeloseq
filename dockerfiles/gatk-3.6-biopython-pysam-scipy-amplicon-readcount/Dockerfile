FROM ubuntu:xenial
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL \
    description="Image for tools used in the CLE"

RUN apt-get update -y && apt-get install -y \
    ant \
    apt-utils \
    default-jdk \
    default-jre \
    git \
    libncurses5-dev \
    ncurses-dev \
    unzip \
    wget \
    zlib1g-dev

RUN apt-get update -y && apt-get install -y python python-pip python-dev build-essential nodejs
RUN pip install --upgrade pip
RUN pip install biopython
RUN pip install pysam
RUN pip install scipy
COPY addAmpliconInfoAndCountReads.py /usr/bin/addAmpliconInfoAndCountReads.py

RUN ln -s /usr/bin/unzip /bin/unzip

#GATK 3.6#
ENV maven_package_name apache-maven-3.3.9
ENV gatk_dir_name gatk-protected
ENV gatk_version 3.6
RUN cd /tmp/ && wget -q http://mirror.nohup.it/apache/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.zip

# LSF: Comment out the oracle.jrockit.jfr.StringConstantPool.
RUN cd /tmp/ \
    && git clone --recursive https://github.com/broadgsa/gatk-protected.git \
    && cd /tmp/gatk-protected && git checkout tags/${gatk_version} \
    && sed -i 's/^import oracle.jrockit.jfr.StringConstantPool;/\/\/import oracle.jrockit.jfr.StringConstantPool;/' ./public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/varianteval/VariantEval.java \
    && mv /tmp/gatk-protected /opt/${gatk_dir_name}-${gatk_version}
RUN cd /opt/ && unzip /tmp/${maven_package_name}-bin.zip \
    && rm -rf /tmp/${maven_package_name}-bin.zip LICENSE NOTICE README.txt \
    && cd /opt/ \
    && cd /opt/${gatk_dir_name}-${gatk_version} && /opt/${maven_package_name}/bin/mvn verify -P\!queue \
    && mv /opt/${gatk_dir_name}-${gatk_version}/protected/gatk-package-distribution/target/gatk-package-distribution-${gatk_version}.jar /opt/GenomeAnalysisTK.jar \
    && rm -rf /opt/${gatk_dir_name}-${gatk_version} /opt/${maven_package_name}

#LDAP and NSS
RUN apt-get update -y && apt-get install -y libnss-sss tzdata
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime

#LSF: Java bug that need to change the /etc/timezone.
#     The above /etc/localtime is not enough.
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

