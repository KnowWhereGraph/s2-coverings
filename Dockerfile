FROM python:3.12

# install some required packages/libraries
RUN apt update && \
    apt install -y --no-install-recommends \
    git \
    make \
    cmake \
    g++ \
    openssl \
    swig \
    python3-dev \
    sed \
    bash \
    gdal-bin \
    libgdal-dev

# Install Abseil at version 20240722.0
RUN cd /tmp && \
    git clone https://github.com/abseil/abseil-cpp.git && \
    cd abseil-cpp && \
    git checkout 20240722.0 && \
    mkdir build && cd build && \
    cmake -DCMAKE_CXX_STANDARD=17 -DCMAKE_POSITION_INDEPENDENT_CODE=ON .. && \
    cmake --build . -j $(nproc) --target install && \
    cd /tmp && \
    rm -rf abseil-cpp

# Install S2 Geometry v0.11.1
RUN cd /tmp && \
    git clone https://github.com/google/s2geometry.git && \
    cd s2geometry && \
    git checkout v0.11.1 && \
    mkdir build &&  cd build && \
    cmake -DCMAKE_CXX_STANDARD=17 -DBUILD_TESTS=OFF -DWITH_PYTHON=ON ..

# Install the SWIG Python bindings
RUN cd /tmp/s2geometry && \
    sed -i '23i "-DBUILD_TESTS=OFF", "-DCMAKE_CXX_STANDARD=17",' setup.py && \
    python3 -m pip install cmake_build_extension wheel && \
    python3 setup.py bdist_wheel && \
    python3 -m pip install "$(find dist -name *.whl)" && \
    cd /tmp && \
    rm -rf s2geometry

# Set the working directory
WORKDIR /s2

# Copy requirements.txt to the working directory inside the image
COPY requirements.txt .

# Install the dependencies
RUN python3 -m ensurepip --upgrade
RUN python3 -m pip install --upgrade setuptools
RUN python3 -m pip install --no-cache-dir -r requirements.txt