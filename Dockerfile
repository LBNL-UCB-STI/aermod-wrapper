# Use the latest Ubuntu image
FROM ubuntu:22.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    make \
    wget \
    m4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    zlib1g-dev \
    vim \
    less \
    curl \
    libnuma1 \
    libncursesw5 \
    libtinfo5 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy the NetCDF fix script into the image
COPY netcdf_3.6_fix.sh /app/netcdf_3.6_fix.sh
RUN chmod +x /app/*.sh

# Install pgf90
COPY nvhpc-21-11_21.11_amd64.deb ./nvhpc-21-11_21.11_amd64.deb
COPY nvhpc-2021_21.11_amd64.deb ./nvhpc-2021_21.11_amd64.deb
RUN apt-get install ./nvhpc-21-11_21.11_amd64.deb ./nvhpc-2021_21.11_amd64.deb
RUN rm ./nvhpc-21-11_21.11_amd64.deb ./nvhpc-2021_21.11_amd64.deb
ENV PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/bin:$PATH

# Install NetCDF 3.6.1
RUN wget https://www.gfd-dennou.org/arch/ucar/netcdf/old/netcdf-3.6.1.tar.gz && \
    tar -xzvf netcdf-3.6.1.tar.gz && \
    bash /app/netcdf_3.6_fix.sh netcdf-3.6.1/src && \
    cd netcdf-3.6.1/src && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf netcdf-3.6.1 netcdf-3.6.1.tar.gz

# Clean up the environment
RUN apt-get purge -y build-essential wget m4 libcurl4-openssl-dev libxml2-dev zlib1g-dev && \
    apt-get autoremove -y && \
    apt-get clean

# Install SSH server
RUN apt-get update && apt-get install -y openssh-server
RUN mkdir /var/run/sshd

# Configure SSH (optional but recommended)
RUN echo 'root:root' | chpasswd
RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

EXPOSE 22

# Copy the AERMOD sources and build script
COPY aermod /app/aermod
COPY aermod_build.sh /app/aermod/build.sh
RUN chmod +x /app/aermod/build.sh

# Compile AERMOD (fortran .f files)
WORKDIR /app/aermod
RUN ./build.sh

# Copy the MMIF sources
COPY mmif /app/mmif

# Compile MMIF (fortran .f90 files)
WORKDIR /app/mmif
RUN make

# Scripts to organize execution
COPY entrypoint.sh /app/scripts/entrypoint.sh
RUN chmod +x /app/scripts/*.sh

# Add apps to PATH
ENV PATH="/app/aermod:/app/mmif:${PATH}"

# Create symlinks for convenience
RUN ln -sf /app/aermod/aermod /usr/local/bin/aermod
RUN ln -sf /app/mmif/mmif /usr/local/bin/mmif

# Set the working directory
WORKDIR /app

# Set entrypoints
ENTRYPOINT ["/app/scripts/entrypoint.sh"]
