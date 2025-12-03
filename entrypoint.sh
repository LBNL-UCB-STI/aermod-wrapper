#!/bin/bash

command_to_execute="${1,,}"
shift;  # will remove first arg from the "$@"
function_args="$@";  # will use all args except first one


# Check the first parameter
case "$command_to_execute" in
    "aermod")
	cd /data
        aermod "$function_args"
        ;;

    "mmif")
	cd /data
        mmif "$function_args"
        ;;

    "ssh") /usr/sbin/sshd -D ;;

    *)
        echo "Usage: aermod|mmif|ssh [arg1 [arg2 [arg3]]]"
        echo "Data folder should be mounted into /data, the application will be executed inside."
        ;;
esac

