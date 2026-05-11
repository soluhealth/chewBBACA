#!/bin/bash
set -euo pipefail

IMAGE="chewbbaca"

while getopts :v:: opt; do
  case ${opt} in
    v)
      RELEASE_VERSION=${OPTARG}
      ;;
    ?)
      echo "Invalid option: -${OPTARG}."
      exit 1
      ;;
  esac
done

VERSIONED_IMAGE=${IMAGE}:${RELEASE_VERSION}

read -p "Release as $VERSIONED_IMAGE (y/n)?" CONT
if [ "$CONT" != "y" ]; then
  echo 'Release canceled'
  exit 1
fi

LOCAL_TAG="solu/$VERSIONED_IMAGE"
echo "Building image locally with tag $LOCAL_TAG"
docker build . --tag $LOCAL_TAG --platform linux/amd64

REMOTE_TAG="us-central1-docker.pkg.dev/solu-platform/solu/$VERSIONED_IMAGE"
echo "Tagging image with $REMOTE_TAG"
docker tag $LOCAL_TAG $REMOTE_TAG

echo "Pushing image to $REMOTE_TAG"
docker push $REMOTE_TAG

echo "Done"
