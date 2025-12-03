#!/bin/bash

ssh_port=22222
aermod_data=/home/local_user/data/aermod
mmif_data=/home/local_user/data/mmif4

docker run --rm -d --privileged \
        -p $ssh_port:22 \
        -v $aermod_data:/data/aermod \
        -v $mmif_data:/data/mmif \
        --name aermod_mmif \
        aermod_mmif:latest \
        ssh

