#!/bin/bash

install_docker_ce() {
    echo "Installing Docker CE"
    sudo apt-get update

    # Allow Apt to work over https
    sudo apt-get install apt-transport-https ca-certificates curl software-properties-common

    # Add docker GPG key
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

    # Add docker repo
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

    # Update the package list and install docker ce
    sudo apt-get update && sudo apt-get install docker-ce -y
}

install_docker_compose() {
    echo "Installing Docker Compose"
    # Get latest release
    sudo curl -L https://github.com/docker/compose/releases/download/1.21.2/docker-compose-$(uname -s)-$(uname -m) -o \
    /usr/local/bin/docker-compose

    # Apply executable permissions to the binary
    sudo chmod +x /usr/local/bin/docker-compose
}

install_google_cloud_sdk() {
    echo "Installing Google Cloud SDK"
    # Create environment variable for correct distribution
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"

    # Add the Cloud SDK distribution URI as a package source
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

    # Import the Google Cloud Platform public key
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -

    # Update the package list and install the Cloud SDK
    sudo apt-get update && sudo apt-get install google-cloud-sdk -y

}

setup_bravo_env() {
    mkdir -p /data/cache/igv_cache
    mkdir -p /data/cram
    mkdir -p /data/genomes
    mkdir -p /data/coverage
    mkdir -p /data/import_vcf
}

install_docker_ce
install_docker_compose
install_google_cloud_sdk
setup_bravo_env
