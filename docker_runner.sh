#!/bin/bash

# Docker runner script for R analysis environment
# This script runs the zjardyn/pelagic_amplicon Docker container
# with the project directory mounted and radian as the default shell

# Set the project directory (current directory)
PROJECT_DIR="/Users/zjardynhood/Desktop/pelagic-microplastics-amplicon"

# Docker image name
IMAGE_NAME="zjardyn/pelagic_amplicon:latest"

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker is not running. Please start Docker first."
    exit 1
fi

# Check if the image exists locally
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo "Image $IMAGE_NAME not found locally. Pulling from Docker Hub..."
    docker pull "$IMAGE_NAME"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to pull image $IMAGE_NAME"
        exit 1
    fi
fi

# Run the container with radian
echo "Starting R environment with radian..."
echo "Project directory: $PROJECT_DIR"
echo "Image: $IMAGE_NAME"
echo ""

docker run -it --rm \
    -v "$PROJECT_DIR:/work" \
    -v "$PROJECT_DIR:$PROJECT_DIR" \
    -w "$PROJECT_DIR" \
    "$IMAGE_NAME" \
    radian
