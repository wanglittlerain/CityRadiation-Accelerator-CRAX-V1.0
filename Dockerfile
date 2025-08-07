FROM gcc:13-bookworm AS build

# Set Arrow version as environment variable
ENV ARROW_VERSION=12.0.1

RUN apt-get update && \
    # add libarrow-dev source
    apt-get install -y -V ca-certificates lsb-release wget && \
    wget https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb && \
    apt-get install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb && \
    apt-get update && \
    apt-get install -y \
    # install build dependencies
    build-essential \
    libssl-dev \
    tar \
    # install crax dependencies
    libarrow-dev=${ARROW_VERSION}* \
    libeigen3-dev

# Create symbolic link for Eigen to match expected path
RUN ln -s /usr/include/eigen3/Eigen /usr/include/Eigen

# Fetch Boost 1.87.0 source
RUN wget https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz && \
    tar xzf boost_1_87_0.tar.gz && \
    mv boost_1_87_0/boost /usr/include/ && \
    rm -rf boost_1_87_0 boost_1_87_0.tar.gz

# Install CMake 3.27.2 with architecture detection
RUN ARCH=$(uname -m) && \
    case ${ARCH} in \
        x86_64) CMAKE_ARCH="linux-x86_64" ;; \
        aarch64|arm64) CMAKE_ARCH="linux-aarch64" ;; \
        *) echo "Unsupported architecture: ${ARCH}, building from source" && \
           wget https://github.com/Kitware/CMake/releases/download/v3.27.2/cmake-3.27.2.tar.gz && \
           tar -xzf cmake-3.27.2.tar.gz && \
           cd cmake-3.27.2 && \
           ./bootstrap && \
           make -j$(nproc) && \
           make install && \
           cd .. && \
           rm -rf cmake-3.27.2 cmake-3.27.2.tar.gz && \
           exit 0 ;; \
    esac && \
    wget https://github.com/Kitware/CMake/releases/download/v3.27.2/cmake-3.27.2-${CMAKE_ARCH}.sh -q -O /tmp/cmake-install.sh && \
    chmod u+x /tmp/cmake-install.sh && \
    mkdir -p /opt/cmake && \
    /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake && \
    ln -s /opt/cmake/bin/* /usr/local/bin/ && \
    rm /tmp/cmake-install.sh

COPY . /CRAX
WORKDIR /CRAX

# only build required binaries
RUN sh build.sh 1

FROM debian:12-slim AS run
COPY --from=build /CRAX/bin /CRAX

FROM scratch AS output
LABEL org.opencontainers.image.source=https://github.com/wanglittlerain/CityRadiation-Accelerator-CRAX
COPY --from=build /CRAX/bin /