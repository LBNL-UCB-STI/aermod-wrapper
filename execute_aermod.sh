#!/bin/bash

echo "Will run AERMOD inside '/data' with following arguments: $@"

cp /app/aermod/aermod /data/aermod
cd /data

chmod +x aermod

echo "Executing AERMOD ..."

./aermod "$@"

