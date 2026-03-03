docker stop carl
docker rm carl

docker buildx build --no-cache . -t carl -f carl_dyncoup.docker
#  && \  docker run --name carl -p 8880:8000 carl

docker exec -it carl /bin/bash 