### Docker Installation Guide

#### Prerequisites

* Administrator/root privileges on your machine
* 64-bit processor with virtualization support (for Windows/macOS)
* At least 4GB RAM (8GB+ recommended)
* Internet connection

#### Installation links 

Windows - [Install Docker Desktop on Windows](https://docs.docker.com/desktop/setup/install/windows-install)

MAC OS - [Install Docker Desktop on Mac](https://docs.docker.com/desktop/setup/install/mac-install/)

Linux - [Install Docker Desktop on Linux](https://docs.docker.com/desktop/setup/install/linux/)


#### Verifying the docker installation

To verify installation of docker run following commands:

```
docker pull hello-world:latest
docker run --rm hello-world:latest
```

The output should contain ‘Hello from Docker!’ text and a short explanation how it works.


### 


### AERMOD

The software and data are available at [Air Quality Dispersion Modeling - Preferred and Recommended Models | US EPA](https://www.epa.gov/scram/air-quality-dispersion-modeling-preferred-and-recommended-models)

MMIF

The software, data available at [Air Quality Dispersion Modeling - Related Model Support Programs | US EPA](https://www.epa.gov/scram/air-quality-dispersion-modeling-related-model-support-programs#mmif)


#### Quick Start Guide

Pull the image

```
docker pull haitamlaarabi/aermod_mmif:1.0
```

Start Container

```
docker run --rm \
  -v <path/to/your/aermod or mmif/data>:/data \
  --name aermod_mmif \
  haitamlaarabi/aermod_mmif:1.0 <aermod or mmif> <config from data folder>
```

Notes:

* `<path/to/your/aermod or mmif/data>` needs to be replaced with a real path, the path connects your local data folder to the container
* Linux/Mac example: `/home/username/projects/data`
* Windows example: `C:\Users\Username\Projects\Data`

Your input files go in the mounted folders, results appear in the same locations.


#### Example

* Download `aermod_test.zip`
* Unzip `aermod_test.zip` to user home folder (`/home/username/aermod_test`)
* Execute:
```
docker run --rm \
	-v /home/username/aermod_test:/data \
    --name aermod_mmif \
    haitamlaarabi/aermod_mmif:1.0 aermod olm.inp
```
* Execution should produce lines like `+Now Processing Data For Day No.  X of 1999`


#### How to build the image

In order to build the docker image: 
* Download nvhpc-2021_21.11_amd64.deb from the [NVIDIA website] (https://developer.nvidia.com/nvidia-hpc-sdk-2111-downloads)
* Download nvhpc-21-11_21.11_amd64.deb from the [NVIDIA website] (https://developer.nvidia.com/nvidia-hpc-sdk-2111-downloads)
* Build docker image: 
```
docker build . -t <image tag>
```

Notes:
- The current docker image contains Aermod v23132. 


